
import gzip
import os
import urllib2
import pysam
import pysamstats
import numpy


#version of numpy used = 1.16.2
#version of pysam used = 0.15.2
#version of pysamstats used = 1.1.2

#Note: this code requires:

# - an indexed reference sequence for Mulatta alpha globin, AlphaGlobRefMulatta1_simpler.fa  (see supplementary material for
#details of the sequence used as a reference sequence).

# - a list of the SRR identifiers of the reads for each animal to be used from the SRA (see supplementary material for details of
#how these were obtained).

# - a list of additional identifiers to generate the correct URLs for each animal We obtained these by trial and error.

# - both the SRR identifiers and additional numerical identifiers are provided in recordedvals.txt.

#Note also lines 60 to 70 have been commented out - these lines need to be uncommented in order to download the
#reads, but since this step only needs to be carried out once, it is likely to be convenient to comment it out again
#before making any alterations to the mapping steps or other elements of the code.


working_folder = '/home/bridget/Desktop/MulattaProject/data/'

os.chdir(working_folder)

beg = '136647'
end = '137050'

fastafilepath = 'All_SupportedSequences_UPDATED.fa'
depthfilepath='AllSupSeqReadDepths.csv'

output_handle = open(fastafilepath, "w")

output_handle2 = open(depthfilepath, 'w')

with open('recordedvals.txt', 'r') as filehandle:
    for line in filehandle:

        SRR_values = line.split(' ')

        SRR_code = SRR_values[0].replace("\r\n", "")

        val1 = SRR_values[1]

        val2 = SRR_values[2].replace("\n", "")

        print SRR_code, val1, val2

        filename = SRR_code + '_' + beg + '_' + end + '.fq.gz'

        # urladdress = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/?path=%2Fnetmnt%2Ftraces04%2Fsra' + val1 + '%2FSRR%2F00' + val2 + '%2F' + SRR_code + '&run=' + SRR_code + '&acc=NC_007877.1&ref=Chr20&range=' + beg + '+-' + end + '&src=0&output=fastq&output_to=File'
        #
        # urlinfo = urllib2.urlopen(urladdress)
        #
        # fa = open(filename, 'w')
        #
        # info = urlinfo.read()
        #
        # fa.write(info)
        #
        # fa.close()

        filenameunzipped_1 = SRR_code + '_' + beg + '_' + end + '_1.fq'

        outputfilehandle1 = open(filenameunzipped_1, 'w')


        with gzip.open(filename) as filehandle1:
            for line1 in filehandle1:
                outputfilehandle1.write(line1)

            outputfilehandle1.close()


        SRR_code_aligned_sorted = SRR_code + '_' + beg + '_' + end + '_aligned_sorted.bam'

        import subprocess

        bowtiestring = 'bowtie2 -x AlphaGlobRefMulatta1_simpler_index -U ' + filenameunzipped_1 + ' | samtools view -bS | samtools sort -o ' + SRR_code_aligned_sorted

        # print bowtiestring

        command1 = [bowtiestring]
        subprocess.call(command1, shell=True, cwd=working_folder)

        indexstring = 'samtools index ' + SRR_code_aligned_sorted
        command2 = [indexstring]
        subprocess.call(command2, shell=True, cwd=working_folder)

        seq_length=204

        Sequence = ['N'] * (seq_length)

        zeros = 0 * (seq_length)

        ReadDepth=['0'] * (seq_length)

        mybam = pysam.AlignmentFile(SRR_code_aligned_sorted)

        StoredVals = numpy.zeros(((seq_length+200), 4))

        CountOverall = 0

        SeqCalled=0

        # iterate over statistics, one record at a time
        for rec in pysamstats.stat_variation(mybam, 'AlphaGlobRefMulatta1_simpler.fa', min_baseq = 20):
            A_count = rec['A']
            T_count = rec['T']
            C_count = rec['C']
            G_count = rec['G']

            StoredVals[rec['pos']][0] = A_count
            StoredVals[rec['pos']][1] = T_count
            StoredVals[rec['pos']][2] = C_count
            StoredVals[rec['pos']][3] = G_count

        for x in range(seq_length):

            A_count = StoredVals[x+100][0]
            T_count = StoredVals[x+100][1]
            C_count = StoredVals[x+100][2]
            G_count = StoredVals[x+100][3]

            TotalCount = A_count + T_count + C_count + G_count

            ReadDepth[x]=TotalCount

            CountOverall = CountOverall + TotalCount

            if TotalCount <= 20 and TotalCount>=4:

                SeqCalled=SeqCalled+1

                if A_count >= 3 and T_count <= 1 and C_count <= 1 and G_count <= 1:
                    Sequence[x] = 'A'

                if C_count >= 3 and A_count <= 1 and T_count <= 1 and G_count <= 1:
                    Sequence[x] = 'C'

                if T_count >= 3 and A_count <= 1 and C_count <= 1 and G_count <= 1:
                    Sequence[x] = 'T'

                if G_count >= 3 and T_count <= 1 and C_count <= 1 and A_count <= 1:
                    Sequence[x] = 'G'

                if A_count >= 2 and T_count >= 2 and C_count <= 1 and G_count <= 1:
                    Sequence[x] = 'W'

                if A_count >= 2 and G_count >= 2 and C_count <= 1 and T_count <= 1:
                    Sequence[x] = 'R'

                if A_count >= 2 and C_count >= 2 and G_count <= 1 and T_count <= 1:
                    Sequence[x] = 'M'

                if T_count >= 2 and G_count >= 2 and C_count <= 1 and A_count <= 1:
                    Sequence[x] = 'K'

                if T_count >= 2 and C_count >= 2 and A_count <= 1 and G_count <= 1:
                    Sequence[x] = 'Y'

                if C_count >= 2 and G_count >= 2 and A_count <= 1 and T_count <= 1:
                    Sequence[x] = 'S'

                if C_count >= 2 and G_count >= 2 and A_count >= 2:
                    Sequence[x] = 'Z'

                if T_count >= 2 and G_count >= 2 and A_count >= 2:
                    Sequence[x] = 'Z'

                if C_count >= 2 and G_count >= 2 and T_count >= 2:
                    Sequence[x] = 'Z'

                if C_count >= 2 and A_count >= 2 and T_count >= 2:
                    Sequence[x] = 'Z'

            if TotalCount > 20:

                SeqCalled = SeqCalled + 1

                if A_count > 0.85*TotalCount and T_count <= 0.05*TotalCount and C_count <= 0.05*TotalCount and G_count <= 0.05*TotalCount:
                    Sequence[x] = 'A'

                if C_count > 0.85*TotalCount and A_count <= 0.05*TotalCount and T_count <= 0.05*TotalCount and G_count <= 0.05*TotalCount:
                    Sequence[x] = 'C'

                if T_count > 0.85*TotalCount and A_count <= 0.05*TotalCount and C_count <= 0.05*TotalCount and G_count <= 0.05*TotalCount:
                    Sequence[x] = 'T'

                if G_count > 0.85*TotalCount and T_count <= 0.05*TotalCount and C_count <= 0.05*TotalCount and A_count <= 0.05*TotalCount:
                    Sequence[x] = 'G'

                if A_count > 0.05*TotalCount and T_count > 0.05*TotalCount and C_count <= 0.05*TotalCount and G_count <= 0.05*TotalCount:
                    Sequence[x] = 'W'

                if A_count > 0.05*TotalCount and G_count > 0.05*TotalCount and C_count <= 0.05*TotalCount and T_count <= 0.05*TotalCount:
                    Sequence[x] = 'R'

                if A_count > 0.05*TotalCount and C_count > 0.05*TotalCount and G_count <= 0.05*TotalCount and T_count <= 0.05*TotalCount:
                    Sequence[x] = 'M'

                if T_count > 0.05*TotalCount and G_count > 0.05*TotalCount and C_count <= 0.05*TotalCount and A_count <= 0.05*TotalCount:
                    Sequence[x] = 'K'

                if T_count > 0.05*TotalCount and C_count > 0.05*TotalCount and A_count <= 0.05*TotalCount and G_count <= 0.05*TotalCount:
                    Sequence[x] = 'Y'

                if C_count > 0.05*TotalCount and G_count > 0.05*TotalCount and A_count <= 0.05*TotalCount and T_count <= 0.05*TotalCount:
                    Sequence[x] = 'S'

                if C_count > 0.05*TotalCount and G_count > 0.05*TotalCount and A_count > 0.05*TotalCount:
                    Sequence[x] = 'Z'

                if T_count > 0.05*TotalCount and G_count > 0.05*TotalCount and A_count > 0.05*TotalCount:
                    Sequence[x] = 'Z'

                if C_count > 0.05*TotalCount and G_count > 0.05*TotalCount and T_count > 0.05*TotalCount:
                    Sequence[x] = 'Z'

                if C_count > 0.05*TotalCount and A_count > 0.05*TotalCount and T_count > 0.05*TotalCount:
                    Sequence[x] = 'Z'



        if SeqCalled > 0:

            SequenceString = ''

            ReadDepthString=''

            for x in range(seq_length):
                SequenceString = SequenceString + Sequence[x]
                ReadDepthString=ReadDepthString+','+ str(ReadDepth[x])

            output_handle.write(">" + SRR_code + '\n' + SequenceString + '\n')

            output_handle2.write(">" + SRR_code + ReadDepthString + '\n')

    output_handle.close()
    output_handle2.close()

