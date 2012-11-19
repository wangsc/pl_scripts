/usr/bin/perl /usr/local/share/perl/5.10.1/ExtUtils/xsubpp  -typemap /usr/share/perl/5.10/ExtUtils/typemap   FOO_07fc.xs > FOO_07fc.xsc && mv FOO_07fc.xsc FOO_07fc.c
cc -c  -I/home/swang/workspace/EA_Lab/test -D_REENTRANT -D_GNU_SOURCE -DDEBIAN -fno-strict-aliasing -pipe -fstack-protector -I/usr/local/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -g   -DVERSION=\"0.00\" -DXS_VERSION=\"0.00\" -fPIC "-I/usr/lib/perl/5.10/CORE"   FOO_07fc.c
FOO_07fc.xs:7:10: warning: missing terminating " character
FOO_07fc.xs: In function ‘greet’:
FOO_07fc.xs:7: error: missing terminating " character
FOO_07fc.xs:8:1: warning: missing terminating " character
FOO_07fc.xs:8: error: missing terminating " character
FOO_07fc.xs:9: error: expected expression before ‘}’ token
FOO_07fc.xs:9: error: expected ‘;’ before ‘}’ token
make: *** [FOO_07fc.o] Error 1
