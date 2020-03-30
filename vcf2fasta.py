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


def filter_vcf_list(vcf_list, min_qual, use_filter_column, min_depth=0, min_maping_qual=0, max_depth=10000):
    filtered_vcf_list = []

    if use_filter_column:
        print('Already filtered')
        for line in vcf_list:
            if line.filter == 'PASS':
                filtered_vcf_list.append(line)

    else:
        for line in vcf_list:
            if line.qual != '.':
                if float(line.qual) > min_qual:
                    filtered_vcf_list.append(line)

    if len(filtered_vcf_list) == 0 and len(vcf_list) > 0:
        print('All variants removed by filter')

    else:
        return filtered_vcf_list


def export_to_fasta_aln(vcf_list_of_objects, out_file_name, ignore_mv_sites, alt_only):
    # print first info line
    number_of_variants = len(vcf_list_of_objects)
    list_of_samples = vcf_list_of_objects[0].sample_dict.keys()
    number_of_samples = len(list_of_samples)

    out_fasta = open(out_file_name + '.fa', 'w')

    if alt_only is False:
        # Export the reference sequence
        out_fasta.write('>' + vcf_list_of_objects[0].chr + '\n')
        nuc_string = ''

        for ref_variant in vcf_list_of_objects:
            nuc_string += ref_variant.refAllele
        out_fasta.write(nuc_string + '\n')

    # Export the other sequences
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
                    if int(variant.sample_dict[a_sample]['AD'].split(',')[0]) < int(
                            variant.sample_dict[a_sample]['AD'].split(',')[1]):

                        if variant.sample_dict[a_sample]['GT'] == './.':
                            nuc_string += variant.refAllele
                        else:
                            print(variant.sample_dict[a_sample])
                            print(variant.altAllele)
                            print(variant.pos)
                        # print(genotype)
                        # print(variant.altAllele.split(',')[genotype - 1])
                        nuc_string += variant.altAllele

                elif len(variant.altAllele) == 1:

                    if ignore_het:
                        if variant.sample_dict[a_sample]['GT'] == '1' or variant.sample_dict[a_sample]['GT'] == '1/1':
                            nuc_string += variant.altAllele
                        if variant.sample_dict[a_sample]['GT'] == '.' or variant.sample_dict[a_sample]['GT'] == '0/1':
                            nuc_string += variant.refAllele

                    else:
                        if variant.sample_dict[a_sample]['GT'] == '1' or \
                                variant.sample_dict[a_sample]['GT'] == '1/1' or \
                                variant.sample_dict[a_sample]['GT'] == '0/1' or \
                                variant.sample_dict[a_sample]['DP'] != '.':
                            nuc_string += variant.altAllele

                        elif variant.sample_dict[a_sample]['GT'] == '.' or variant.sample_dict[a_sample]['GT'] == './.':
                            # Bug here, what is the point of this?
                            nuc_string += variant.refAllele
                        else:
                            print(variant.sample_dict[a_sample]['GT'])

        out_fasta.write(nuc_string + '\n')


def main(inFile, output_dir, qual_filter, ignore_mv_sites, use_filter_column, ignore_het, alt_only):
    """

    :param inFile:
    :param output_dir:
    :return:
    """
    invcf = input_parser(inFile)

    filtered_vcf_list = filter_vcf_list(invcf, qual_filter, use_filter_column)

    export_to_fasta_aln(filtered_vcf_list, output_dir, ignore_mv_sites, alt_only)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''A tool for the conversion of vcf files into aligned fasta format so that they may be used in 
        phylogenetic tree construction.
        ''', epilog="""Version 0.1""")

    parser.add_argument('-i', '--input_file', type=str, help='Input file in vcf format')
    parser.add_argument('-o', '--output_file', type=str, help='Output fasta file name')
    parser.add_argument('-q', '--quality', type=int, default=2000, help='Quality threshold for variants')
    parser.add_argument('-m', '--ignore_mv_sites', type=bool, default=True,
                        help='Exclude sites with more than one alt allele')
    parser.add_argument('-f', '--already_filtered', type=bool, default=False,
                        help='Use the FILTER column to keep / drop variants')
    parser.add_argument('-e', '--ignore_heterozygous', type=bool, default=False,
                        help='Set to True to ignore variants that are 0/1')
    parser.add_argument('-a', '--alt_only', type=bool, default=False,
                        help='Set to True to only include the alt alleles in the output fasta')

    args = parser.parse_args()

    inFile = args.input_file
    output_dir = args.output_file
    qual_filter = args.quality
    ignore_mv_sites = args.ignore_mv_sites
    use_filter_column = args.already_filtered
    ignore_het = args.ignore_heterozygous
    alt_only = args.alt_only

    main(inFile, output_dir, qual_filter, ignore_mv_sites, use_filter_column, ignore_het, alt_only)
