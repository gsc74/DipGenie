/* Narrow library adapter for pinned HTSbox pileup.c, MIT (see
 * third_party/HTSBOX-LICENSE.txt). Upstream source remains unmodified.
 * The C++ caller serializes getopt access. Stream state is thread-local;
 * process stdout and random-number state are never changed. */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "sam.h"
static _Thread_local FILE *dg_pileup_output;
static _Thread_local int dg_pileup_read_error;
static int dg_checked_bam_read(BGZF *fp,bam1_t *b) {
    const int r=bam_read1(fp,b);
    if(r < -1)dg_pileup_read_error=1;
    return r;
}
#undef stdout
#define stdout dg_pileup_output
#define printf(...) fprintf(dg_pileup_output,__VA_ARGS__)
#define putchar(c) fputc((c),dg_pileup_output)
#define puts(s) (fputs((s),dg_pileup_output)<0 ? EOF : fputc('\n',dg_pileup_output))
// Only the deterministic VCF path is exposed, never the random-consensus mode.
#define srand48(seed) ((void)(seed))
#define bam_read1 dg_checked_bam_read
#define main_pileup dg_upstream_pileup
#include "../third_party/htsbox/pileup.c"
#undef main_pileup
#undef bam_read1
#undef srand48
#undef puts
#undef putchar
#undef printf
#undef stdout

int dg_embedded_pileup(int argc,char **argv,FILE *output) {
    int result,old_optind=optind,old_opterr=opterr,old_optopt=optopt;
    char *old_optarg=optarg;
    dg_pileup_output=output;dg_pileup_read_error=0;
    optind=0;opterr=0;
    result=dg_upstream_pileup(argc,argv);
    if(dg_pileup_read_error||fflush(output)||ferror(output))result=-1;
    dg_pileup_output=NULL;
    optind=old_optind;opterr=old_opterr;optopt=old_optopt;optarg=old_optarg;
    return result;
}
