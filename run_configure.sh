./configure --prefix=/home1/st18003/scratch/software/install --disable-gsl \
CPPFLAGS="-I /home1/st18003/scratch/projects/a4md/_install/include \
-I/software/python/3.6.3/include/python3.6m -I/software/python/3.6.3/include/python3.6m  -Wno-unused-result -Wsign-compare  -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes" \
LDFLAGS="-L /home1/st18003/scratch/projects/a4md/_install/lib -lingest -lretrieve -la4md_cmn -ldspaces -ldscommon -ldart -lm \
-L /software/boost/1.68-gcc-4.8.5/lib -lboost_serialization \
-L/software/python/3.6.3/lib -lpython3.6m -lpthread -ldl  -lutil  -Xlinker -export-dynamic"
