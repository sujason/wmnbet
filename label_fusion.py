#!/usr/bin/env python
"""
Generic label fusion script
JSON-format for configuration:
{
    'path': 'data/',
    'images': 'image.name.nii.gz',
    'labels': {
        'brain': ['', ''],
    }
    'template': 'template.nii.gz',
    'postprocess': 'c3d %s %s',

}
If template doesn't exist, register all priors
Postprocess is a command to apply to labels afterwards, %s should exist twice
"""
import os
import sys
import glob
import time
import json
import argparse
import tempfile
import parallel
from shutil import rmtree
from functools import partial
from collections import OrderedDict
from datetime import timedelta
from libraries.imgtools import check_run, check_warps, sanitize_input, flip_lr, label_average, label_fusion_majority, label_fusion_picsl, ants_compose_a_to_b, ants_apply_only_warp
from libraries.ants_nonlinear import ants_nonlinear_registration, bias_correct


parser = argparse.ArgumentParser(description='Perform label fusion given a data path of priors and configuration JSON file.')
subparsers = parser.add_subparsers(help='sub-command help', dest='action')
check_parser = subparsers.add_parser('check', help='check json') 
check_parser.add_argument('config_json', help='configuration JSON file')
main_parser = subparsers.add_parser('do', help='do registration and label fusion')
main_parser.add_argument('config_json', help='configuration JSON file')
main_parser.add_argument('input_image', help='input NIfTI image')
main_parser.add_argument('output_path', help='the output file for single ROI or directory for multiple ROIs')
main_parser.add_argument('roi_names', metavar='roi_names', nargs='+', help='a space separated list of one or more ROIs')
main_parser.add_argument('-w', '--warp', metavar='path', help='looks for {path}InverseWarp.nii.gz and {path}Affine.txt instead of basing it off input_image.')
main_parser.add_argument('-F', '--forcereg', action='store_true', help='force ANTS registration to WMnMPRAGE mean brain template. The --warp argument can be then used to specify the output path.')
main_parser.add_argument('-p', '--processes', nargs='?', default=None, const=None, type=int, help='number of parallel processes to use.  If unspecified, automatically set to number of CPUs.')
main_parser.add_argument('-v', '--verbose', action='store_true', help='verbose mode')
main_parser.add_argument('-R', '--right', action='store_true', help='segment right thalamus')
main_parser.add_argument('--labelfusion', choices=('majority', 'average', 'correlation', 'jointfusion', 'STEPS'), default='majority', help='label fusion method')
main_parser.add_argument('--tempdir', help='temporary directory to store registered atlases.  This will not be deleted as usual.')


exec_options = {'echo': False, 'suppress': True}


def do_nothing(*args):
    pass


def populate_priors(config_json, verbose=False):
    with open(config_json, 'r') as f:
        config = json.load(f)
    error_flag = False
    if not os.path.isabs(config['path']):
        config['path'] = os.path.join(os.path.normpath(os.path.dirname(config_json)), config['path'])
    prior_subjects = filter(os.path.isdir, glob.glob(os.path.join(config['path'], '*')))
    priors = {}
    for p in prior_subjects:
        name = os.path.basename(p)
        priors[name] = {}
        if verbose:
            print name
        for label_name, candidates in config['labels'].iteritems():
            priors[name].setdefault('labels', {})[label_name] = filter(os.path.exists, [os.path.join(p, el) for el in candidates])
            if verbose:
                print '\t', label_name, len(priors[name]['labels'][label_name]), priors[name]['labels'][label_name]
        find_files = ['image']
        if 'template' in config:
            find_files += ['template_transform_affine', 'template_transform_warp', 'template_transform_inversewarp']
        for el in find_files:
            priors[name][el] = os.path.join(p, config[el])
            if not os.path.exists(priors[name][el]):
                print '%s is missing %s' % (name, el)
                error_flag = True
    if 'template' in config:
        template = os.path.join(config['path'], config['template'])
        if not os.path.exists(template):
            print 'Template is missing: %s' % template
            error_flag = True
    else:
        template = None
    if not error_flag:
        print 'All good!'
    else:
        print 'Problems were found!'
    return config, template, priors, error_flag


def warp_priors_via_template(priors, labels, input_image, input_transform_prefix, output_path, pool, execute):
    """
    Warp labels from a set of priors to the input_image via an intermediate template space.
    """
    # Combine warps
    cmds = []
    combined_warps = {}
    for subject, datum in priors.iteritems():
        # Name ['image', 'labels', 'template_transform_inversewarp', 'template_transform_affine', 'template_transform_warp']
        subj_output_path = os.path.join(output_path, subject)
        try:
            os.mkdir(subj_output_path)
        except OSError:
            # Exists
            pass
        a_transform_prefix = datum['template_transform_affine'].replace('Affine.txt', '')
        combined_warps[subject] = combined_warp = os.path.join(subj_output_path, 'Warp.nii.gz')
        if not os.path.exists(combined_warp):
            _, cmd = ants_compose_a_to_b(
                a_transform_prefix,
                b_path=input_image,
                b_transform_prefix=input_transform_prefix,
                output=combined_warp,
                execute=do_nothing,
            )
            cmds.append(cmd)
    pool.map(execute, cmds)
    # Warp images and labels
    cmds = []
    output_labels = {}
    for subject, datum in priors.iteritems():
        output_image = os.path.join(output_path, subject, os.path.basename(datum['image']))
        # Anatomical image
        if not os.path.exists(output_image):
            _, cmd = ants_apply_only_warp(
                template=input_image,
                input_image=datum['image'],
                input_warp=combined_warps[subject],
                output_image=output_image,
                switches='--use-BSpline',
                execute=do_nothing,
            )
            cmds.append(cmd)
        # ROIs
        for label in labels:
            label_fnames = datum['labels'][label]
            # We allow multiple candidates for a given ROI per subject
            for label_fname in label_fnames:
                warped_label = os.path.join(output_path, subject, os.path.basename(label_fname))
                output_labels.setdefault(label, []).append((warped_label, output_image))
                if not os.path.exists(warped_label):
                    _, cmd = ants_apply_only_warp(
                        template=input_image,
                        input_image=label_fname,
                        input_warp=combined_warps[subject],
                        output_image=warped_label,
                        switches='--use-NN',
                        execute=do_nothing,
                    )
                    cmds.append(cmd)
    pool.map(execute, cmds)
    return output_labels


def main(args, config, template, priors, temp_path, pool, exec_options):
    parallel_command = partial(parallel.command, **exec_options)
    input_image = orig_input_image = args.input_image
    output_path = args.output_path
    labels = []
    for roi in args.roi_names:
        if roi not in config['labels'].keys():
            print '%s is in valid, must be one of these: %s' % (roi, ', '.join(config['labels'].keys()))
            sys.exit(1)
        labels.append(roi)

    if args.warp:
        warp_path = args.warp
    else:
        head, tail = os.path.split(input_image)
        tail = tail.replace('.nii', '').replace('.gz', '')
        warp_path = os.path.join(temp_path, tail)

    t = time.time()
    # FSL automatically converts .nii to .nii.gz
    sanitized_image = os.path.join(temp_path, os.path.basename(input_image) + ('.gz' if input_image.endswith('.nii') else ''))
    print '--- Reorienting image. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
    if not os.path.exists(sanitized_image):
        input_image = sanitize_input(input_image, sanitized_image, parallel_command)
        if args.right:
            print '--- Flipping along L-R. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
            flip_lr(input_image, input_image, parallel_command)
        print '--- Correcting bias. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
        bias_correct(input_image, input_image, **exec_options)
    else:
        print 'Skipped, using %s' % sanitized_image
        input_image = sanitized_image

    if template:
        print '--- Registering to mean brain template. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
        if args.forcereg or not check_warps(warp_path):
            if args.warp:
                print 'Saving output as %s' % warp_path
            else:
                warp_path = os.path.join(temp_path, tail)
                print 'Saving output to temporary path.'
            ants_nonlinear_registration(template, input_image, warp_path, **exec_options)
        else:
            print 'Skipped, using %sInverseWarp.nii.gz and %sAffine.txt' % (warp_path, warp_path)
        print '--- Warping prior labels and images. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
        # TODO should probably use output from warp_atlas_subject instead of hard coding paths in create_atlas
        # TODO make this more parallel
        warped_labels = warp_priors_via_template(
            priors=priors,
            labels=labels,
            input_image=input_image,
            input_transform_prefix=warp_path,
            output_path=temp_path,
            pool=pool,
            execute=parallel_command,
        )
    else:
        # register each prior
        raise NotImplementedError('registration of each prior directly, instead of template')
    # print '--- Forming subject-registered atlases. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
    # atlases = pool.map(partial(create_atlas, path=temp_path, subjects=subjects, target='', echo=exec_options['echo']),
    #     [{'label': label, 'output_atlas': os.path.join(temp_path, label+'_atlas.nii.gz')} for label in warped_labels])
    # atlases = dict(zip(warped_labels, zip(*atlases)[0]))
    # atlas_image = atlases['WMnMPRAGE_bias_corr']
    # atlas_images = warped_labels['WMnMPRAGE_bias_corr']
    print '--- Performing label fusion. --- Elapsed: %s' % timedelta(seconds=time.time() - t)
    if args.labelfusion == 'majority':
        pool.map(label_fusion_majority,
            [dict(
                atlas_labels=[roi for roi, img in warped_labels[label]],
                output_label=os.path.join(temp_path, label+'.nii.gz'),
                execute=parallel_command
            ) for label in labels])
    elif args.labelfusion == 'average':
        pool.map(label_average,
            [dict(
                atlas_labels=[roi for roi, img in warped_labels[label]],
                output_label=os.path.join(temp_path, label+'.nii.gz'),
                execute=parallel_command
            ) for label in labels])
    else:
        # pool.map(partial(label_fusion_picsl, input_image, atlas_images),
        #     [dict(
        #         atlas_labels=[roi for roi, img in warped_labels[label]],
        #         output_label=os.path.join(temp_path, label+'.nii.gz'),
        #         rp=optimal_picsl[label]['rp'],
        #         rs=optimal_picsl[label]['rs'],
        #         beta=optimal_picsl[label]['beta'],
        #         **exec_options
        #     ) for label in labels])
        raise NotImplementedError(args.labelfusion)
    files = [(os.path.join(temp_path, label + '.nii.gz'), os.path.join(output_path, label + '.nii.gz')) for label in labels]
    if args.right:
        pool.map(flip_lr, files)
        files = [(out_file, out_file) for in_file, out_file in files]
    if 'postprocess' in config:
        print '--- Postprocessing labels. --- Elapsed: %s' % timedelta(seconds=time.time() - t)
        pool.map(parallel_command, [config['postprocess'] % (in_file, out_file) for in_file, out_file in files])
        files = [(out_file, out_file) for in_file, out_file in files]
    # Resort output to original ordering
    pool.map(parallel_command,
        ['%s %s %s %s' % ('/RAID/THOMAS/distro/swapdimlike.py', in_file, orig_input_image, out_file) for in_file, out_file in files])
    print '--- Finished --- Elapsed: %s' % timedelta(seconds=time.time() - t)


if __name__ == '__main__':
    args = parser.parse_args()
    # print args
    if args.action == 'check':
        populate_priors(args.config_json, True)
        sys.exit(0)
    else:
        config, template, priors, error_flag = populate_priors(args.config_json)
        if error_flag:
            print 'Quitting due to missing files.'
            sys.exit(1)
    if args.verbose:
        exec_options['echo'] = True
    pool = parallel.BetterPool(args.processes)
    print 'Running with %d processes.' % pool._processes
    # TODO Add path of script to command()
    # os.environ['PATH'] += os.pathsep + os.path.abspath(os.path.dirname(sys.argv[0]))
    if args.tempdir:
        temp_path = args.tempdir
        if not os.path.exists(temp_path):
            print 'Making %s' % os.path.abspath(temp_path)
            os.makedirs(temp_path)
    else:
        temp_path = tempfile.mkdtemp(dir=os.path.dirname(args.output_path))
    try:
        main(args, config, template, priors, temp_path, pool, exec_options)
    finally:
        if not args.tempdir:
            try:
                rmtree(temp_path)
            except OSError as exc:
                if exc.errno != 2:  # Code 2 - no such file or directory
                    raise
