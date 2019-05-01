from __future__ import division, print_function

def cwd_safe_to_write(do_not_overwrite_files):
  '''Check if the current working directory is safe to write, i.e. doesn't contain files specified by the caller'''
  import glob
  for f1 in do_not_overwrite_files:
    for f2 in glob.glob(f1):
      from libtbx.utils import Sorry
      sorry_msg = "Cannot overwrite existing file: " + f2
      raise Sorry(sorry_msg)

if __name__=="__main__":
  import sys
  # can test one or more comma-separated filenames with wildcards
  if len(sys.argv) > 1:
    cwd_safe_to_write(sys.argv[1].split(','))
  else:
    print ("Usage: safe_to_write <file1>,<file2>,...")

