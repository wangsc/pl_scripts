/usr/bin/perl /usr/local/share/perl/5.10.1/ExtUtils/xsubpp  -typemap /usr/share/perl/5.10/ExtUtils/typemap   FOO_6950.xs > FOO_6950.xsc && mv FOO_6950.xsc FOO_6950.c
cc -c  -I/home/swang/workspace/EA_Lab/test -D_REENTRANT -D_GNU_SOURCE -DDEBIAN -fno-strict-aliasing -pipe -fstack-protector -I/usr/local/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -g   -DVERSION=\"0.00\" -DXS_VERSION=\"0.00\" -fPIC "-I/usr/lib/perl/5.10/CORE"   FOO_6950.c
FOO_6950.xs:7:10: warning: missing terminating " character
FOO_6950.xs: In function ‘greet’:
FOO_6950.xs:7: error: missing terminating " character
FOO_6950.xs:8: error: expected expression before ‘}’ token
FOO_6950.xs:8: error: expected ‘;’ before ‘}’ token
make: *** [FOO_6950.o] Error 1
