import subprocess as sb

#Will take 6.5 hours, mostly from pipeline1 and 2

sb.call('python3 pipeline1.py', shell=True)
sb.call('python3 pipeline2.py', shell=True)
sb.call('python3 make_raw_data_file.py', shell=True)
sb.call('python3 making_annotations.py', shell=True)
sb.call('python3 make_csv_file.py', shell=True)
sb.call('python3 make_tsv_file.py', shell=True)
