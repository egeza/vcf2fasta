

import argparse


class Variant:
	def __init__(self, vcf_line, headder_line):
		vcf_line = vcf_line.strip()
		vcf_line = vcf_line.split('\t')
		self.chr = vcf_line[0]
		self.pos = vcf_line[1]
		self.id = vcf_line[2]
		self.refAllele = vcf_line[3]
		self.altAllele = vcf_line[4]
		self.qual = vcf_line[5]
		self.filter = vcf_line[6]
		self.info = vcf_line[7]
		self.format = vcf_line[8]

		format_index = vcf_line[8].strip()
		format_index = format_index.split(':')

		sample_dict = {}
		count = 9
		for a_sample in headder_line:
			sample_info = dict(zip(format_index, vcf_line[count].split(':')))
			sample_dict[a_sample] = sample_info
			count += 1

		self.sample_dict = sample_dict


def input_parser(file_path):

	if file_path[-4:] == ".vcf":
		list_of_objs = []
		# Deal with random info at start
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if line[0:2] != "##" and line[0] == "#":
				vcf_headder_line = line.split('\t')
				vcf_headder_line[0] = vcf_headder_line[0][1:]
				vcf_headder_line[-1] = vcf_headder_line[-1].strip()
				vcf_samples = vcf_headder_line[9:]

			if not line.startswith('#'):
				variant_record = Variant(line, vcf_samples)

				list_of_objs.append(variant_record)

		return list_of_objs


def filter_vcf_list(vcf_list, min_qual, min_depth=0, min_maping_qual=0, max_depth=10000):
	filtered_vcf_list = []

	for line in vcf_list:
		if line.qual != '.':
			if float(line.qual) > min_qual:
				filtered_vcf_list.append(line)

	return filtered_vcf_list

def export_to_fasta_aln(vcf_list_of_objects, out_file_name):
	# print first info line
	number_of_variants = len(vcf_list_of_objects)
	list_of_samples = vcf_list_of_objects[0].sample_dict.keys()
	number_of_samples = len(list_of_samples)

	out_fasta = open(out_file_name + '.fa', 'w')

	for a_sample in list_of_samples:
		out_fasta.write('>' + a_sample + '\n')
		nuc_string = ''
		for variant in vcf_list_of_objects:
			if len(variant.sample_dict[a_sample].keys()) < 2 and len(variant.altAllele) == 1:
				nuc_string += variant.refAllele
			else:
				# Account for multiple alt alleles found at a site
				if len(variant.altAllele) != 1 and ignore_mv_sites is False:
					print("resolving multi allele sites")
					if int(variant.sample_dict[a_sample]['AD'].split(',')[0]) < int(variant.sample_dict[a_sample]['AD'].split(',')[1]):
						#print(variant.altAllele)
						#print(variant.sample_dict[a_sample])

						#print(variant.sample_dict[a_sample]['AD'].split(','))
						#print(variant.qual)
						#print(variant.sample_dict[a_sample])
						if variant.sample_dict[a_sample]['GT'] == './.':
							nuc_string += variant.refAllele
						else:
							print(variant.sample_dict[a_sample])
							print(variant.altAllele)
							print(variant.pos)
						#print(genotype)
						#print(variant.altAllele.split(',')[genotype - 1])
						nuc_string += variant.altAllele

				elif len(variant.altAllele) == 1:
					nuc_string += variant.altAllele

		out_fasta.write(nuc_string + '\n')


def main(inFile, output_dir, qual_filter=2000):
	"""

	:param inFile:
	:param output_dir:
	:return:
	"""
	ignore_mv_sites = True

	invcf = input_parser(inFile)

	filtered_vcf_list = filter_vcf_list(invcf, qual_filter)

	export_to_fasta_aln(filtered_vcf_list, output_dir)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='''A de-multiplexer tool Input for primers is csv separated by a comma. 
		Headings include: gene_name,fwd_sequence,fwd_PID_len,rev_sequence,rev_pid,overlapping
		overlapping = True if the the fwd and rev reads overlap, else overlapping = False
		''',  epilog="""Version 0.1""")

	parser.add_argument('-i', '--input_file', type=str, help='Input file in vcf format')
	parser.add_argument('-o', '--output_file', type=str, help='Output fasta file name')
	parser.add_argument('-q', '--quality', type=int, default=2000, help='Quality threshold for variants')
	parser.add_argument('-m', '--ignore_mv_sites', type=int, default=2000, help='Exclude sites with more than one alt allele')

	args = parser.parse_args()

	inFile = args.input_file
	output_dir = args.output_file
	main(inFile, output_dir)

