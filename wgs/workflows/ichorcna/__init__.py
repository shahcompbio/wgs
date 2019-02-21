import os
import tasks
import pypeliner
import pypeliner.managed as mgd
from pypeliner_utils import helpers


def create_hmmcopy_workflow(bams, out_dir, config):
	workflow = pypeliner.workflow.Workflow()

	out_dir = os.path.join(out_dir, '{sample_id}')

	workflow.setobj(
		obj=mgd.OutputChunks('sample_id'),
		value=bams.keys(),
	)


	workflow.transform(
		name='hmmcopy_readcounter',
		func=tasks.hmmcopy_readcounter,
		axes=('sample_id',),
		ctx={
			'mem':	config['memory']['low'],
			'pool_id':	config['pools']['standard'],
			'ncpus': 1,
		},
		args=(
			mgd.InputFile('input.bam', 'sample_id', fnames=bams),
			mgd.TempOutputFile('infile.wig', 'sample_id'),
			config,
		)
	)


	workflow.transform(
		name='run_ichorcna',
		func=tasks.run_hmmcopy,
		axes=('sample_id',),
		ctx={
			'mem':	config['memory']['med'],
			'pool_id':	config['pools']['standard'],
			'ncpus': 1,
		},
		args=(
			mgd.TempInputFile('infile_copy.obj', 'sample_id'),
			mgd.TempInputFile('infile_copy.txt', 'sample_id'),
			mgd.TempOutputFile('hmmcopy_res.obj', 'sample_id'),
			mgd.TempOutputFile('hmmcopy_segments.txt', 'sample_id'),
			mgd.OutputFile(tumour_table_out, 'sample_id'),
			mgd.InputInstance('sample_id'),
			config,
		)
	)


    return workflow