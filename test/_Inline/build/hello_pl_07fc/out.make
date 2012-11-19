/usr/bin/perl /usr/local/share/perl/5.10.1/ExtUtils/xsubpp  -typemap /usr/share/perl/5.10/ExtUtils/typemap   hello_pl_07fc.xs > hello_pl_07fc.xsc && mv hello_pl_07fc.xsc hello_pl_07fc.c
cc -c  -I/home/swang/workspace/EA_Lab/test -D_REENTRANT -D_GNU_SOURCE -DDEBIAN -fno-strict-aliasing -pipe -fstack-protector -I/usr/local/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -g   -DVERSION=\"0.00\" -DXS_VERSION=\"0.00\" -fPIC "-I/usr/lib/perl/5.10/CORE"   hello_pl_07fc.c
hello_pl_07fc.xs:7:10: warning: missing terminating " character
hello_pl_07fc.xs: In function ‘greet’:
hello_pl_07fc.xs:7: error: missing terminating " character
hello_pl_07fc.xs:8:1: warning: missing terminating " character
hello_pl_07fc.xs:8: error: missing terminating " character
hello_pl_07fc.xs:9: error: expected expression before ‘}’ token
hello_pl_07fc.xs:9: error: expected ‘;’ before ‘}’ token
make: *** [hello_pl_07fc.o] Error 1
