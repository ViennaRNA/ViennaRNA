#include <stdio.h>
#include <g2.h>

#ifdef DO_PS
#include <g2_PS.h>
#endif
#ifdef DO_FIG
#include <g2_FIG.h>
#endif
#ifdef DO_X11
#include <g2_X11.h>
#endif
#ifdef DO_GD
#include <g2_gd.h>
#endif
#ifdef DO_WIN32
#include <g2_win32.h>
#endif
#ifdef DO_WMF32
#include <g2_win32.h>
#endif

#define maxdev 10

int main()
{
    int d, dev[maxdev]={-1, -1, -1, -1, -1,-1, -1, -1, -1, -1};
    int ndev=0;

    printf("\nG2_VERSION %s\n", G2_VERSION);
    
    d=g2_open_vd();      /* open virtual device */

    printf("Adding..");
    
#ifdef DO_PS
    printf("..PS");
    dev[ndev]=g2_open_PS("g2_arc.ps", g2_A4, g2_PS_land);
    g2_attach(d, dev[ndev]);
    ndev++;
 
    printf("..EPSF");
    dev[ndev]=g2_open_EPSF("g2_arc.eps");
    g2_attach(d, dev[ndev]);
    ndev++;
    printf("..EPSF_CLIP");
    dev[ndev]=g2_open_EPSF_CLIP("g2_arc_clip.eps",400,400);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_FIG
    printf("..FIG");
    dev[ndev]=g2_open_FIG("g2_arc.fig");
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_X11
    printf("..X11");
    dev[ndev]=g2_open_X11(400,400);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_GD
    printf("..GD(png)");
    dev[ndev]=g2_open_gd("g2_arc.png", 400, 400,g2_gd_png);
    g2_attach(d, dev[ndev]);
    ndev++;
    printf("..GD(jpeg)");
    dev[ndev]=g2_open_gd("g2_arc.jpg", 400, 400,g2_gd_jpeg);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_WIN32
    printf("..WIN32");
    dev[ndev]=g2_open_win32(400, 400,"g2_arc",0);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_WMF32
    printf("..WMF32");
    dev[ndev]=g2_open_win32(775, 575,"g2_arc.emf",1);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
    g2_set_auto_flush(d,0);

    printf("\n");

    g2_string(d,10, 50,"90,90");
    g2_string(d,10,100,"90,360");
    g2_string(d,10,150,"90,0");
    g2_string(d,10,200,"360,90");
    g2_string(d,10,250,"-45,45");
    g2_string(d,10,300,"-45,-315");
    g2_string(d,10,350,"-495,-405");

    g2_arc(d, 100,  50, 45, 25, 90, 90);
    g2_arc(d, 100, 100, 45, 25, 90, 360);
    g2_arc(d, 100, 150, 45, 25, 90, 0);
    g2_arc(d, 100, 200, 45, 25, 360, 90);
    g2_arc(d, 100, 250, 45, 25, -45, 45);
    g2_arc(d, 100, 300, 45, 25, -45, -315);
    g2_arc(d, 100, 350, 45, 25, -495, -405);

    g2_filled_arc(d, 100,  50, 40, 20, 90, 90);
    g2_filled_arc(d, 100, 100, 40, 20, 90, 360);
    g2_filled_arc(d, 100, 150, 40, 20, 90, 0);
    g2_filled_arc(d, 100, 200, 40, 20, 360, 90);
    g2_filled_arc(d, 100, 250, 40, 20, -45, 45);
    g2_filled_arc(d, 100, 300, 40, 20, -45, -315);
    g2_filled_arc(d, 100, 350, 40, 20, -495, -405);


    
    g2_string(d,320, 50,"0,0");
    g2_string(d,320,100,"360,90");
    g2_string(d,320,150,"0,90");
    g2_string(d,320,200,"90,360");
    g2_string(d,320,250,"45,-45");
    g2_string(d,320,300,"-315,-45");
    g2_string(d,320,350,"-405,-495");

    g2_arc(d, 250,  50, 45, 25, 0, 0);
    g2_arc(d, 250, 100, 45, 25, 360, 90);
    g2_arc(d, 250, 150, 45, 25, 0, 90);
    g2_arc(d, 250, 200, 45, 25, 90, 360);
    g2_arc(d, 250, 250, 45, 25, 45, -45);
    g2_arc(d, 250, 300, 45, 25, -315, -45);
    g2_arc(d, 250, 350, 45, 25, -405, -495);

    g2_filled_arc(d, 250,  50, 40, 20, 90, 90);
    g2_filled_arc(d, 250, 100, 40, 20, 360, 90);
    g2_filled_arc(d, 250, 150, 40, 20, 0, 90);
    g2_filled_arc(d, 250, 200, 40, 20, 90, 360);
    g2_filled_arc(d, 250, 250, 40, 20, 45, -45);
    g2_filled_arc(d, 250, 300, 40, 20, -315, -45);
    g2_filled_arc(d, 250, 350, 40, 20, -405, -495);

    g2_flush(d);
    printf("\nDone.\n[Enter]\n");
    getchar();
    g2_close(d);
    return 0;
}

