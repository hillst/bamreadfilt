# Bug seems to originate from reads which overlap multiple variants being included twice in teh bam.
#   Suggestions:
#       adding another hashset that looks for read_ids (assuming the read pair has a unique id) or having a custom writer that handles it.
#       creating a seperate routine that works after the bam has been writen to remove duplicates 
#       creating a separate program that works after the bam has been written to remove duplicates
#       somehow make sure the bam pointer only moves forward ? (coordinate sorted?)

BAM=test/test.bam
VCF=test/test4k.vcf

./bamreadfilt --bam ${BAM} --vcf ${VCF}
samtools sort new_out.bam > new_out.sorted.bam
samtools index new_out.sorted.bam

#run again on the output of the last step
./bamreadfilt --bam new_out.sorted.bam --vcf ${VCF}
samtools sort new_out.bam > new_out.again.bam
samtools index new_out.again.bam

samtools view new_out.again.bam > new_out.again.sam
samtools view new_out.sorted.bam > new_out.sorted.sam

#they are different at this point, key difference being the EXACT SAME READ EVERYTHING occurs multiple times.
diff new_out.again.sam new_out.sorted.sam

## how do we resolve this? we either need to track the reads and apply a uniqueness filter or we need to write an 
##      an extension to the write function that enforces uniqueness (with an internal set of somesort)
