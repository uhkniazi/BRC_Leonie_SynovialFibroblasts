# File: header.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 17/10/2017
# Desc: header file with paths and global settings

gcsWD = getwd()
g_did = 9

f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}