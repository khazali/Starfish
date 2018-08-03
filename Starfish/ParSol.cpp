#include "Globals.h"
#ifndef PETSC_ENABLED

FloatType ParaSolver(FloatType *A, int *row, int *col, FloatType *b, FloatType *xx) {
	LocalVector<double> x;
	LocalVector<double> rhs;
	LocalMatrix<double> mat;


	mat.AllocateCSR("A", TotalMSize, PoreNO, PoreNO);
	x.Allocate("x", PoreNO);
	rhs.Allocate("rhs", PoreNO);

	mat.CopyFromCSR(row, col, A);
	x.CopyFromData(xx);
	rhs.CopyFromData(b);

	/*mat.SetDataPtrCSR(&prow, &pcol, &pA, "A", TotalMSize, PoreNO, PoreNO);
	x.SetDataPtr(&pxx, "x", PoreNO);
	rhs.SetDataPtr(&pb, "rhs", PoreNO);*/

	mat.MoveToAccelerator();
	x.MoveToAccelerator();
	rhs.MoveToAccelerator();

	// Linear Solver
	CG<LocalMatrix<double>, LocalVector<double>, double > ls;
	TNS<LocalMatrix<double >, LocalVector<double >, double > p;

	ls.Init(1e-50, 1e-20, 1e100, 50000);


	ls.SetOperator(mat);
	ls.SetPreconditioner(p);

	ls.Build();

	ls.Verbose(0);

	//mat.info();	

	ls.Solve(rhs, &x);
	x.CopyToData(xx);
	ls.Clear();
	mat.Clear();
	x.Clear();
	rhs.Clear();

	return 0;
}
#endif
