The barcode file is a tab delimited unix '\n' line ended file. It has one line per library/sample, with a header line.
The file consists of 6 columns, the columns are as follows:
1. Library/sample name, this is the name of the resultant file for this barcode, a .fq will be automatically appended. Please do not use special characters like '-','.',',' or spaces. Also, try not to use the same name twice.

2. The read type. This program can handle most combinations of barcoded, single indexed and dual indexed. The read type is as follows:
Read Type:
0 - Single ended, barcoded, no index
1 - Pair ended, barcoded, no index
2 - Single ended, one index
3 - Pair ended, one index
4 - Pair ended, two indexes
5 - Single ended, two index.
Please note that while the default operating mode for indexed read is to use a separate index file, it will look at the end of the read name line after the '#' character if no index file(s) are specified.

3. The first barcode/index

4. The second barcode/index if present, if not use a '.' as a placeholder.

5. The overhang after the barcode. Some barcode schemes have a 'T' or other overhang after the barcode, indicate the overhang here if present.

6. The removed bases from start of reads. IT is sometimes necessary to remove a number of bases from the start of the read after barcode/index matching. The number of bases would be specified here.

For an example of this type of file, please see the sample-barcode-file.txt.
Please use library and file naming and formatting requirements or data will be
lost.


