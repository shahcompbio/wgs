from pypeliner_utils import vcfutils

output_vcf = {
	'1': '/genesis/shahlab/pwalters/MAS-13/titan-temp/tmp/sample_id/T8/titan/chr1.vcf.tmp',
	'2': '/genesis/shahlab/pwalters/MAS-13/titan-temp/tmp/sample_id/T8/titan/chr2.vcf.tmp',
	'3': '/genesis/shahlab/pwalters/MAS-13/titan-temp/tmp/sample_id/T8/titan/chr3.vcf.tmp',
	'4': '/genesis/shahlab/pwalters/MAS-13/titan-temp/tmp/sample_id/T8/titan/chr4.vcf.tmp',
	'5': '/genesis/shahlab/pwalters/MAS-13/titan-temp/tmp/sample_id/T8/titan/chr5.vcf.tmp',
}

out = '/genesis/shahlab/pwalters/MAS-13/museq.vcf'

vcfutils.concatenate_vcf(output_vcf, out)