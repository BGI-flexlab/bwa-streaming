cmake_minimum_required(VERSION 3.6)
project(bwa)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "-g -Wall -Wno-unused-function -O2 -DHAVE_PTHREAD -DUSE_MALLOC_WRAPPERS")

set(LIBBWA_SOURCE_FILES
        utils.c
        kthread.c
        kstring.c
        ksw.c
        bwt.c
        bntseq.c
        bwa.c
        bwamem.c
        bwamem_pair.c
        bwamem_extra.c
        malloc_wrap.c
        QSufSort.c
        bwt_gen.c
        rope.c
        rle.c
        is.c
        bwtindex.c)

set(MAIN_SOURCE_FILES
        bamlite.c
        bwa.c
        bwape.c
        bwase.c
        bwaseqio.c
        bwashm.c
        bwt_lite.c
        bwtaln.c
        bwtalnpe.c
        bwtalnse.c
        bwtgap.c
        bwtindex.c
        bwtsw2_aux.c
        bwtsw2_chain.c
        bwtsw2_core.c
        bwtsw2_main.c
        bwtsw2_pair.c
        fastmap.c
        kopen.c
        main.c
        maxk.c
        pemerge.c  filter.c filter_main.c)

add_library(bwa ${LIBBWA_SOURCE_FILES})
add_executable(bwa-streaming ${MAIN_SOURCE_FILES})
target_link_libraries(bwa-streaming m z pthread bwa)