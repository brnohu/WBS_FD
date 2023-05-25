///***************************************************************************
// File:        FRAMEWRK.H
//
// Purpose:     Header file for Framework library
//
// Platform:    Microsoft Windows
//
// Comments:
//              Include this file in any source files
//              that use the framework library.
//
// From the Microsoft Excel Developer's Kit, Version 14
// Copyright (c) 1997 - 2010 Microsoft Corporation. All rights reserved.
///***************************************************************************



#define rwMaxO8			(65536)
#define colMaxO8		(256)
#define cchMaxStz		(255)
#define MAXSHORTINT		0x7fff

// 
// Function prototypes
//

#ifdef __cplusplus
extern "C" {
#endif

	void far cdecl debugPrintf(LPSTR lpFormat, ...);
	LPSTR GetTempMemory(size_t cBytes);													// Savine
	void FreeAllTempMemory(void);														// Savine
	//int cdecl Excel(int xlfn, LPXLOPER pxResult, int count, ...);
	LPXLOPER TempNum(double d);
	LPXLOPER TempStr(LPSTR lpstr);
	LPXLOPER TempBool(int b);
	LPXLOPER TempInt(short int i);
	//LPXLOPER TempActiveRef(WORD rwFirst,WORD rwLast,BYTE colFirst,BYTE colLast);
	//LPXLOPER TempActiveCell(WORD rw, BYTE col);
	//LPXLOPER TempActiveRow(WORD rw);
	//LPXLOPER TempActiveColumn(BYTE col);
	LPXLOPER TempErr(WORD i);
	LPXLOPER TempMissing(void);
	void InitFramework(void);
	void QuitFramework(void);

	LPXLOPER TempStrConst(const LPSTR lpstr);
	LPXLOPER12 TempNum12(double d);														// Savine
	LPXLOPER12 TempStr12(const XCHAR* lpstr);											// Savine
	LPXLOPER12 TempBool12(BOOL b);
	LPXLOPER12 TempInt12(int i);
	LPXLOPER12 TempActiveRef12(RW rwFirst,RW rwLast,COL colFirst,COL colLast);
	LPXLOPER12 TempActiveCell12(RW rw, COL col);
	LPXLOPER12 TempActiveRow12(RW rw);
	LPXLOPER12 TempActiveColumn12(COL col);
	LPXLOPER12 TempErr12(int i);														// Savine
	LPXLOPER12 TempMissing12(void);
	int cdecl Excel12f(int xlfn, LPXLOPER12 pxResult, int count, ...);					// Savine

	void FreeXLOperT(LPXLOPER pxloper);
	void FreeXLOper12T(LPXLOPER12 pxloper12);
	BOOL ConvertXLRefToXLRef12(LPXLREF pxref, LPXLREF12 pxref12);
	BOOL ConvertXLRef12ToXLRef(LPXLREF12 pxref12, LPXLREF pxref);
	BOOL XLOperToXLOper12(LPXLOPER pxloper, LPXLOPER12 pxloper12);
	BOOL XLOper12ToXLOper(LPXLOPER12 pxloper12, LPXLOPER pxloper);

#ifdef __cplusplus
}
#endif 
