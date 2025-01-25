#ifndef _POLARFLY_TABLES_HPP_
#define _POLARFLY_TABLES_HPP_

#define USE_TABLE_7x3
//#define USE_TABLE_13x4 
//#define USE_TABLE_57x8

#ifdef USE_TABLE_7x3
const int POLARFLY_TABLE_ROWS = 7;
const int POLARFLY_TABLE_COLS = 3;
const int polarfly_connection_table[7][3] = {
    {3,4,0},
    {5,4,1},
    {6,4,2},
    {6,5,0},
    {0,1,2},
    {6,3,1},
    {5,3,2}
};

#elif defined(USE_TABLE_13x4)
const int POLARFLY_TABLE_ROWS = 13;
const int POLARFLY_TABLE_COLS = 4;
const int polarfly_connection_table[13][4] = {
    {4,5,6,0},
    {7,8,6,1},
    {7,5,9,2},
    {4,8,9,3},
    {7,10,0,3},
    {11,8,0,2},
    {12,9,0,1},
    {4,10,1,2},
    {11,5,1,3},
    {12,6,2,3},
    {12,11,7,4},
    {12,10,8,5},
    {11,10,9,6}
};

#elif defined(USE_TABLE_57x8)
const int POLARFLY_TABLE_ROWS = 57;
const int POLARFLY_TABLE_COLS = 8;
const int polarfly_connection_table[57][8] = {
    {8,9,10,11,12,13,14,0},
    {15,16,17,11,18,19,20,1},
    {21,22,23,24,12,19,25,2},
    {26,27,28,24,18,13,29,3},
    {26,22,17,30,31,32,14,4},
    {21,27,10,33,34,32,20,5},
    {15,9,23,33,31,35,29,6},
    {8,16,28,30,34,35,25,7},
    {26,36,23,0,37,38,7,20},
    {39,22,40,0,34,18,6,41},
    {42,43,28,0,31,5,19,44},
    {45,24,33,46,0,1,47,30},
    {48,16,49,0,2,50,32,29},
    {21,51,17,0,3,35,52,53},
    {15,27,54,0,55,4,56,25},
    {21,36,28,1,55,50,6,14},
    {39,27,49,1,31,12,7,53},
    {48,51,23,1,34,4,13,44},
    {42,9,40,1,3,38,32,25},
    {26,43,10,1,2,35,56,41},
    {8,22,54,1,37,5,52,29},
    {15,36,40,30,2,5,13,53},
    {39,9,28,46,2,4,52,20},
    {8,27,17,47,2,38,6,44},
    {45,11,37,2,34,31,3,55},
    {42,51,54,33,2,18,7,14},
    {8,36,49,33,3,4,19,41},
    {39,16,23,47,3,5,56,14},
    {15,22,10,46,3,50,7,44},
    {48,43,54,30,3,12,6,20},
    {21,43,40,47,11,4,7,29},
    {42,16,10,24,37,4,6,53},
    {45,35,50,4,18,12,5,38},
    {26,51,49,46,11,5,6,25},
    {48,9,17,24,55,5,7,41},
    {45,32,13,52,6,7,56,19},
    {45,39,42,21,15,8,26,48},
    {8,51,40,24,31,50,56,20},
    {8,43,23,46,55,18,32,53},
    {45,36,43,9,22,27,16,51},
    {21,9,49,30,37,18,56,44},
    {26,9,54,47,34,50,19,53},
    {48,36,10,47,31,18,52,25},
    {39,51,10,30,55,38,19,29},
    {45,54,17,28,40,49,23,10},
    {39,36,54,24,11,35,32,44},
    {48,22,28,33,11,38,56,53},
    {42,27,23,30,11,50,52,41},
    {42,36,17,46,34,12,56,29},
    {26,16,40,33,55,12,52,44},
    {15,51,28,47,37,12,32,41},
    {39,43,17,33,37,50,13,25},
    {42,22,49,47,55,35,13,20},
    {21,16,54,46,31,38,13,41},
    {45,44,25,20,41,53,14,29},
    {15,43,49,24,34,38,52,14},
    {48,27,40,46,37,35,19,14}
};

#else
#error "No table size selected! Uncomment one of the USE_TABLE_* defines"
#endif

#endif // _POLARFLY_TABLES_HPP_
