#ifndef highOrderData_H
#define highOrderData_H

#include "labelList.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// 结构化网格模板数据 (在文件作用域声明, 由 structuredStencilBuild.H 初始化)
// 这些变量在 hy2Foam_include.H 中通过 #include "highOrderData.H" 引入,
// 使它们在整个 hy2Foam.C 的函数作用域内可见

// 轴对齐面方向列表: 0=非轴对齐, 1=+x, -1=-x, 2=+y, -2=-y, 3=+z, -3=-z
extern Foam::labelList axisFaces;

// 轴对齐面个数
extern Foam::label nAxisFaces;

// 每个单元的 6 方向邻居, -1 表示边界
extern Foam::labelList cell_plus_x;
extern Foam::labelList cell_minus_x;
extern Foam::labelList cell_plus_y;
extern Foam::labelList cell_minus_y;
extern Foam::labelList cell_plus_z;
extern Foam::labelList cell_minus_z;

// 是否结构化网格
extern bool isStructuredGrid;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
