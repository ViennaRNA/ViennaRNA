/*------------------------------------------------------------------------------
#
#
#	MultiFinder-Aware Simple Input/Output Window resource
#
#	for ReadSeq
#
------------------------------------------------------------------------------*/

#include "systypes.r"
#include "types.r"


resource 'MENU' (20000, preload) {
	20000,
	textMenuProc,
	0x7FFFFFFD,
	enabled,
	apple,
	{	/* array: 2 elements */
		/* [1] */
		"About ReadSeq…", noIcon, noKey, noMark, plain,
		/* [2] */
		"-", noIcon, noKey, noMark, plain
	}
};

resource 'MENU' (20001, preload) {
	20001,
	textMenuProc,
	0x0,
	enabled,
	"File",
	{	/* array: 11 elements */
		/* [1] */
		"New", noIcon, "N", noMark, plain,
		/* [2] */
		"Open", noIcon, "O", noMark, plain,
		/* [3] */
		"-", noIcon, noKey, noMark, plain,
		/* [4] */
		"Close", noIcon, "W", noMark, plain,
		/* [5] */
		"Save", noIcon, "S", noMark, plain,
		/* [6] */
		"Save As…", noIcon, noKey, noMark, plain,
		/* [7] */
		"-", noIcon, noKey, noMark, plain,
		/* [8] */
		"Page Setup…", noIcon, noKey, noMark, plain,
		/* [9] */
		"Print…", noIcon, noKey, noMark, plain,
		/* [10] */
		"-", noIcon, noKey, noMark, plain,
		/* [11] */
		"Quit", noIcon, "Q", noMark, plain
	}
};

resource 'MENU' (20002, preload) {
	20002,
	textMenuProc,
	0x0,
	enabled,
	"Edit",
	{	/* array: 6 elements */
		/* [1] */
		"Undo", noIcon, "Z", noMark, plain,
		/* [2] */
		"-", noIcon, noKey, noMark, plain,
		/* [3] */
		"Cut", noIcon, "X", noMark, plain,
		/* [4] */
		"Copy", noIcon, "C", noMark, plain,
		/* [5] */
		"Paste", noIcon, "V", noMark, plain,
		/* [6] */
		"Clear", noIcon, noKey, noMark, plain
	}
};

resource 'MENU' (20003, preload) {
	20003,
	textMenuProc,
	allEnabled,
	enabled,
	"Font",
	{	/* array: 0 elements */
	}
};

resource 'ALRT' (20000, purgeable) {
	{98, 108, 314, 405},
	20000,
	{	/* array: 4 elements */
		/* [1] */
		OK, visible, silent,
		/* [2] */
		OK, visible, silent,
		/* [3] */
		OK, visible, silent,
		/* [4] */
		OK, visible, silent
	}
};

resource 'ALRT' (20001, purgeable) {
	{40, 20, 150, 260},
	20001,
	{	/* array: 4 elements */
		/* [1] */
		OK, visible, silent,
		/* [2] */
		OK, visible, silent,
		/* [3] */
		OK, visible, silent,
		/* [4] */
		OK, visible, silent
	}
};

resource 'ALRT' (20002, preload) {
	{72, 64, 212, 372},
	20002,
	{	/* array: 4 elements */
		/* [1] */
		OK, visible, silent,
		/* [2] */
		OK, visible, silent,
		/* [3] */
		OK, visible, silent,
		/* [4] */
		OK, visible, silent
	}
};

resource 'DITL' (20000, purgeable) {
	{	/* array DITLarray: 8 elements */
		/* [1] */
		{191, 98, 211, 178},
		Button {
			enabled,
			"OK"
		},
		/* [2] */
		{110, 24, 130, 256},
		StaticText {
			disabled,
			" Copyright © 1990 by d.g.gilbert\n"
		},
		/* [3] */
		{6, 93, 24, 281},
		StaticText {
			disabled,
			"A tool for molecular biology."
		},
		/* [4] */
		{31, 25, 86, 281},
		StaticText {
			disabled,
			"Reads and writes nucleic or protein sequ"
			"ences in various formats. Data files may"
			" have multiple sequences."
		},
		/* [5] */
		{6, 17, 22, 92},
		StaticText {
			disabled,
			"ReadSeq"
		},
		/* [6] */
		{150, 28, 186, 262},
		StaticText {
			disabled,
			"land mail: biology dept., indiana univer"
			"sity, bloomington, in 47405\n"
		},
		/* [7] */
		{129, 25, 153, 258},
		StaticText {
			disabled,
			" e-mail: gilbertd@bio.indiana.edu\n"
		},
		/* [8] */
		{86, 12, 107, 281},
		StaticText {
			disabled,
			"This program may be freely distributed."
		}
	}
};

resource 'DITL' (20001, purgeable) {
	{	/* array DITLarray: 3 elements */
		/* [1] */
		{80, 150, 100, 230},
		Button {
			enabled,
			"OK"
		},
		/* [2] */
		{10, 60, 60, 230},
		StaticText {
			disabled,
			"Error. ^0."
		},
		/* [3] */
		{8, 8, 40, 40},
		Icon {
			disabled,
			2
		}
	}
};

resource 'DITL' (20002, preload) {
	{	/* array DITLarray: 4 elements */
		/* [1] */
		{58, 25, 76, 99},
		Button {
			enabled,
			"Yes"
		},
		/* [2] */
		{86, 25, 104, 99},
		Button {
			enabled,
			"No"
		},
		/* [3] */
		{12, 20, 45, 277},
		StaticText {
			disabled,
			"Save changes before closing?"
		},
		/* [4] */
		{86, 195, 104, 269},
		Button {
			enabled,
			"Cancel"
		}
	}
};

resource 'CNTL' (20000, purgeable, preload) {
	{-1, 465, 272, 481},
	0,
	invisible,
	0,
	0,
	scrollBarProc,
	0,
	""
};

resource 'CNTL' (20001, purgeable, preload) {
	{271, -1, 287, 466},
	0,
	invisible,
	0,
	0,
	scrollBarProc,
	0,
	""
};

data 'pzza' (128, purgeable) {
	$"4D50 5320"                                          /* MPS  */
};

resource 'MBAR' (20000, preload) {
	{	/* array MenuArray: 4 elements */
		/* [1] */
		20000,
		/* [2] */
		20001,
		/* [3] */
		20002,
		/* [4] */
		20003
	}
};

resource 'WIND' (20000, purgeable, preload) {
	{0, 0, 286, 480},
	zoomDocProc,
	invisible,
	noGoAway,
	0x0,
	"untitled"
};

resource 'STR#' (20000, purgeable) {
	{	/* array StringArray: 11 elements */
		/* [1] */
		"You must run on 512Ke or later",
		/* [2] */
		"Application Memory Size is too small",
		/* [3] */
		"Not enough memory to run SIOW",
		/* [4] */
		"Not enough memory to do Cut",
		/* [5] */
		"Cannot do Cut",
		/* [6] */
		"Cannot do Copy",
		/* [7] */
		"Cannot exceed 32,000 characters with Pas"
		"te",
		/* [8] */
		"Not enough memory to do Paste",
		/* [9] */
		"Cannot create window",
		/* [10] */
		"Cannot exceed 32,000 characters",
		/* [11] */
		"Cannot do PasteFont not found"
	}
};

resource 'SIZE' (-1) {
	reserved,
	acceptSuspendResumeEvents,
	reserved,
	canBackground,
	multiFinderAware,
	backgroundAndForeground,
	dontGetFrontClicks,
	ignoreChildDiedEvents,
	not32BitCompatible,
	notHighLevelEventAware,
	onlyLocalHLEvents,
	notStationeryAware,
	dontUseTextEditServices,
	reserved,
	reserved,
	reserved,
	124928,
	38912
};

resource 'SIZE' (0) {
	reserved,
	acceptSuspendResumeEvents,
	reserved,
	canBackground,
	multiFinderAware,
	backgroundAndForeground,
	dontGetFrontClicks,
	ignoreChildDiedEvents,
	not32BitCompatible,
	notHighLevelEventAware,
	onlyLocalHLEvents,
	notStationeryAware,
	dontUseTextEditServices,
	reserved,
	reserved,
	reserved,
	256000,
	38912
};

data 'siow' (0) {
	$"0F52 6561 6453 6571 2069 6E20 5349 4F57"            /* .ReadSeq in SIOW */
};

resource 'BNDL' (128) {
	'siow',
	0,
	{	/* array TypeArray: 2 elements */
		/* [1] */
		'ICN#',
		{	/* array IDArray: 1 elements */
			/* [1] */
			0, 128
		},
		/* [2] */
		'FREF',
		{	/* array IDArray: 1 elements */
			/* [1] */
			0, 128
		}
	}
};

resource 'FREF' (128) {
	'APPL',
	0,
	""
};

resource 'ICN#' (128) {
	{	/* array: 2 elements */
		/* [1] */
		$"0000 0000 0000 0000 0010 4100 0010 2200"
		$"0020 2200 0020 2100 0020 4100 0010 4200"
		$"0010 4200 0010 2200 0020 2100 0020 0100"
		$"00FF FF00 03FF FFE0 0791 03F0 0ED1 0E7C"
		$"1C31 321C 380D C10E 3FFF FFFE 3003 C106"
		$"380D 300E 1E31 0E3C 1FC1 01F8 07FF FFE0"
		$"00FF FE",
		/* [2] */
		$"0000 0000 0000 0000 0010 4100 0010 2200"
		$"0020 2200 0020 2100 0020 4100 0010 4200"
		$"0010 4200 0010 2200 0020 2100 0020 0100"
		$"00FF FF00 03FF FFE0 07FF FFF0 0FFF FFFC"
		$"1FFF FFFC 3FFF FFFE 3FFF FFFE 3FFF FFFE"
		$"3FFF FFFE 1FFF FFFC 1FFF FFF8 07FF FFE0"
		$"00FF FE"
	}
};

