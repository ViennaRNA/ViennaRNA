/*****************************************************************************
**  Copyright (C) 1998-2001  Ljubomir Milanovic & Horst Wagner
**  This file is part of the g2 library
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
******************************************************************************/
#ifndef _G2_PS_DEFINITIONS_H
#define _G2_PS_DEFINITIONS_H

/*
 *
 * Sizes for paper defined in g2_PS.h
 * Size is in 1/72 inch (=0.351mm)
 */
static int g2_PS_paper_size[][2]={
 { 2384, 3370 },   /* g2_A0               -  A0                */
 { 1684, 2384 },   /* g2_A1               -  A1                */
 { 1191, 1684 },   /* g2_A2               -  A2                */
 {  842, 1191 },   /* g2_A3               -  A3                */
 {  595,  842 },   /* g2_A4               -  A4                */
 {  420,  595 },   /* g2_A5               -  A5                */
 {  297,  420 },   /* g2_A6               -  A6                */
 {  210,  297 },   /* g2_A7               -  A7                */
 {  148,  210 },   /* g2_A8               -  A8                */
 {  105,  148 },   /* g2_A9               -  A9                */
 { 2920, 4127 },   /* g2_B0               -  B0                */
 { 2064, 2920 },   /* g2_B1               -  B1                */
 { 1460, 2064 },   /* g2_B2               -  B2                */
 { 1032, 1460 },   /* g2_B3               -  B3                */
 {  729, 1032 },   /* g2_B4               -  B4                */
 {  516,  729 },   /* g2_B5               -  B5                */
 {  363,  516 },   /* g2_B6               -  B6                */
 {  258,  363 },   /* g2_B7               -  B7                */
 {  181,  258 },   /* g2_B8               -  B8                */
 {  127,  181 },   /* g2_B9               -  B9                */
 {   91,  127 },   /* g2_B10              -  B10               */
 {  297,  684 },   /* g2_Comm_10_Envelope -  Comm #10 Envelope */
 {  461,  648 },   /* g2_C5_Envelope      -  C5 Envelope       */
 {  312,  624 },   /* g2_DL_Envelope      -  DL Envelope       */
 {  595,  935 },   /* g2_Folio            -  Folio             */
 {  522,  756 },   /* g2_Executive        -  Executive         */
 {  612,  792 },   /* g2_Letter           -  Letter            */
 {  612, 1008 },   /* g2_Legal            -  Legal             */
 { 1224,  792 },   /* g2_Ledger           -  Ledger            */
 {  792, 1224 }    /* g2_Tabloid          -  Tabloid           */
};


/*
 *
 * PS operators
 *
 */
char *g2_PS_operators[]={
    " /L { lineto } def",			       /* lineto */
    " /St { stroke } def",			       /* stroke */
    " /M { moveto } def",			       /* moveto */
    " /P {",					       /* plot */
    " gsave newpath [] 0 setdash 1 setlinecap 0 setlinewidth",
    " 0.2 sub exch 0.2 sub exch moveto 0.4 0.4 rlineto",
    " stroke grestore} def",
    " /T {",					       /* triangle */
    " newpath",
    " moveto lineto lineto",
    " closepath stroke} def",
    " /FT {",					       /* filled triangle */
    " newpath",
    " moveto lineto lineto",
    " closepath fill} def",
    " /R {",					       /* rectangle */
    " newpath",
    " 3 index 1 index 6 4 roll 5 index 1 index",
    " moveto lineto lineto lineto closepath stroke} def",
    " /FR {",					       /* filled rectangle */
    " newpath",
    " 3 index 1 index 6 4 roll 5 index 1 index",
    " moveto lineto lineto lineto closepath fill} def",
    " /A {",					       /* arc */
    " gsave /g2_old_matrix matrix currentmatrix def newpath",
    " translate scale 0 0 1 5 3 roll arc",
    " g2_old_matrix setmatrix stroke grestore } def",
    " /FA {",					       /* filled arc */
    " gsave /g2_old_matrix matrix currentmatrix def newpath",
    " translate scale 0 0 moveto 0 0 1 5 3 roll arc closepath",
    " g2_old_matrix setmatrix fill grestore } def",
    " /S {",					       /* draw string */
    " gsave newpath",
    " translate 0 0 moveto show",
    " stroke grestore} def",
    "\n",
    NULL
};

#endif /* _G2_PS_DEFINITIONS_H */
