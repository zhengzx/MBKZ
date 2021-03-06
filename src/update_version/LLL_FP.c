
static
long MBKZ(mat_ZZ& BB, mat_ZZ* UU, double delta,
	long beta, long prune, LLLCheckFct check, int ncs, char* qr2)
{
	
	double comp_sum = 0;
	double comp_st1;
	double comp_st2;


	long m = BB.NumRows();
	long n = BB.NumCols();
	long m_orig = m;

	long i, j;
	ZZ MU;

	double t1;
	ZZ T1;
	double *tp;

	init_red_fudge();

	mat_ZZ B;
	B = BB;

	B.SetDims(m + 1, n);


	double **B1;  // approximates B

	typedef double *doubleptr;

	B1 = NTL_NEW_OP doubleptr[m + 2];
	if (!B1) Error("BKZ_FP: out of memory");

	for (i = 1; i <= m + 1; i++) {
		B1[i] = NTL_NEW_OP double[n + 1];
		if (!B1[i]) Error("BKZ_FP: out of memory");
	}

	double **mu;
	mu = NTL_NEW_OP doubleptr[m + 2];
	if (!mu) Error("LLL_FP: out of memory");

	for (i = 1; i <= m + 1; i++) {
		mu[i] = NTL_NEW_OP double[m + 1];
		if (!mu[i]) Error("BKZ_FP: out of memory");
	}


	double *c; // squared lengths of Gramm-Schmidt basis vectors

	c = NTL_NEW_OP double[m + 2];
	if (!c) Error("BKZ_FP: out of memory");

	double *b; // squared lengths of basis vectors

	b = NTL_NEW_OP double[m + 2];
	if (!b) Error("BKZ_FP: out of memory");

	double cbar;

	double *ctilda;
	ctilda = NTL_NEW_OP double[m + 2];
	if (!ctilda) Error("BKZ_FP: out of memory");

	double *vvec;
	vvec = NTL_NEW_OP double[m + 2];
	if (!vvec) Error("BKZ_FP: out of memory");

	double *yvec;
	yvec = NTL_NEW_OP double[m + 2];
	if (!yvec) Error("BKZ_FP: out of memory");

	double *uvec;
	uvec = NTL_NEW_OP double[m + 2];
	if (!uvec) Error("BKZ_FP: out of memory");

	double *utildavec;
	utildavec = NTL_NEW_OP double[m + 2];
	if (!utildavec) Error("BKZ_FP: out of memory");


	long *Deltavec;
	Deltavec = NTL_NEW_OP long[m + 2];
	if (!Deltavec) Error("BKZ_FP: out of memory");

	long *deltavec;
	deltavec = NTL_NEW_OP long[m + 2];
	if (!deltavec) Error("BKZ_FP: out of memory");



	double *sn;//store n vectors'norm during orth_enum
	sn = NTL_NEW_OP double[n + 2];

	double **snv; //store n vectors during orth_enum
	snv = NTL_NEW_OP double*[n + 2];

	for (int i = 1; i <= m + 1; i++)
		snv[i] = NTL_NEW_OP double[n + 2];

	double ***tch; 
	tch = NTL_NEW_OP double**[2];

	for (int i = 0; i <= 1; i++)
		tch[i] = NTL_NEW_OP double*[n + 2];

	for (int j = 0; j <= 1; j++)
		for (int i = 0; i <= m + 1; i++)
			tch[j][i] = NTL_NEW_OP double[n + 2];
	
	int qtch = 0;
	int* stch;
	stch = NTL_NEW_OP int[2];


	//stch, qtch and tch work together to store initial values for orth_enum, only a small part of tch and stch is used in the current version, others are used for multi-thread version. 


	double **un; 
	un = NTL_NEW_OP double*[n + 2];

	for (int i = 1; i <= m + 1; i++)
		un[i] = NTL_NEW_OP double[n + 2];

	double *ylen;
	ylen = NTL_NEW_OP double[m + 2];

	//un and ylen are used to store intermediate results to avoid repeated calculation during the depth first search 

	double **xz;  //store all choices for x_i with  size 5 at most, act as poss_v 
	xz = NTL_NEW_OP double*[5];

	for (int i = 0; i <= 4; i++)
		xz[i] = NTL_NEW_OP double[n + 2];
	

	double *des; //store the average values
	des = NTL_NEW_OP double[n + 2];

	double *jw;//act as poss_v_cnt to store the number of the choices
	jw = NTL_NEW_OP double[n + 2];


	mat_ZZ Ulocal;
	mat_ZZ *U;

	if (UU) {
		Ulocal.SetDims(m + 1, m);
		for (i = 1; i <= m; i++)
			conv(Ulocal(i, i), 1);
		U = &Ulocal;
	}
	else
		U = 0;

	long quit;
	long new_m;
	long z, jj, kk;
	long s, t;
	long h;
	double eta;
	int shuchu = 0;
	char shc[30] = "Basis_";	
	strcpy(shc+6,qr2);	
	char shctra[30] = "trace_";
	strcpy(shctra+6,qr2);
	char shcb[30] = "best_";
	strcpy(shcb+5,qr2);
	char shct[30] = "t";
	strcpy(shct+1,qr2);   

	
	ofstream comp(shct, ios::app);
	ofstream ol;
	ofstream tra1(shctra,ios::app);
	ofstream tra2(shcb);
	tra2.close();	

	comp_st1 = GetTime();

	double ts1;
	double ts2;
	double ts3;
	double ts4;
	double ts5;
	double tmp;
	
	int iii;
	int cs;
	int hl = 0;
	int scsc = 1;
	
	//initialize variables for intermediate result outputs and define temporary variables

	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++) {
			conv(B1[i][j], B(i, j));
			CheckFinite(&B1[i][j]);
		}


	for (i = 1; i <= m; i++) {
		b[i] = InnerProduct(B1[i], B1[i], n);
		CheckFinite(&b[i]);
	}



	m = ll_LLL_FP(B, U, delta, 0, check, B1, mu, b, c, m, 1, quit);
	

	double sr_out;
	int pda;
	double tt;

	double enum_time = 0;
	unsigned long NumIterations = 0;
	unsigned long NumTrivial = 0;
	unsigned long NumNonTrivial = 0;
	unsigned long NumNoOps = 0;
	double tt1;
	long verb = verbose;	
	verbose = 0;

	long clean = 1;

	if (m < m_orig) {
		for (i = m_orig + 1; i >= m + 2; i--) {
			// swap i, i-1

			swap(B(i), B(i - 1));
			if (U) swap((*U)(i), (*U)(i - 1));
		}
	}

	if (!quit && m > 1) {
		if (beta > m) beta = m;

		if (prune > 0)
			ComputeBKZConstant(beta, prune);

		z = 0;
		jj = 0;

		while (z < m - 1) {
			jj++;
			kk = min(jj + beta - 1, m);

			if (jj == m) {
				jj = 1;
				kk = beta;
				clean = 1;
			}

			if (verb) {
				tt = GetTime();
				if (tt > LastTime + LLLStatusInterval)
					BKZStatus(tt, enum_time, NumIterations, NumTrivial,
						NumNonTrivial, NumNoOps, m, B);
			}


			// ENUM

			if (verb) {
				tt1 = GetTime();
			}


			if (prune > 0)
				ComputeBKZThresh(&c[jj], kk - jj + 1);
		
			if (jj == 1) {		
				hl++;
				if(hl%10==0)
				{
					tra2.open(shcb);
					tra2 << "enum: " << jj << " " << kk << " " << c[jj] << " " << powf(c[jj], 0.5);
					tra2 << "   " << hl << endl;
					tra2.close();
				}
				scsc = 1;

				//hl stores how many rounds (run orth_enum once or run tran_enum n times) of enum have been conducted, current best will be output every 10 rounds

			}

			shuchu = 0;
			if (!(hl%2==0&&jj==1)) {// run tran_enum

				cbar = c[jj];
				utildavec[jj] = uvec[jj] = 1;

				yvec[jj] = vvec[jj] = 0;
				Deltavec[jj] = 0;


				s = t = jj;
				deltavec[jj] = 1;

				for (i = jj + 1; i <= kk + 1; i++) {
					ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
					Deltavec[i] = 0;
					vvec[i] = 0;
					deltavec[i] = 1;
				}

				long enum_cnt = 0;

				while (t <= kk) {
					if (verb) {
						enum_cnt++;
						if (enum_cnt > 100000) {
							enum_cnt = 0;
							tt = GetTime();
							if (tt > LastTime + LLLStatusInterval) {
								enum_time += tt - tt1;
								tt1 = tt;
								BKZStatus(tt, enum_time, NumIterations, NumTrivial,
									NumNonTrivial, NumNoOps, m, B);
							}
						}
					}

					ctilda[t] = ctilda[t + 1] +
						(yvec[t] + utildavec[t])*(yvec[t] + utildavec[t])*c[t];

					ForceToMem(&ctilda[t]);  // prevents an infinite loop

					if (prune > 0 && t > jj) {
						eta = BKZThresh(t - jj);
					}
					else
						eta = 0;
					if (ctilda[t] < cbar *((1.0*(beta + jj - t)) / (1.0*(beta)))) { // variable pruning strategy for tran_enum
						//if (ctilda[t] < cbar - eta) {
						if (t > jj) {
							t--;
							t1 = 0;
							for (i = t + 1; i <= s; i++)
								t1 += utildavec[i] * mu[i][t];
							yvec[t] = t1;
							t1 = -t1;
							if (t1 >= 0)
								t1 = ceil(t1 - 0.5);
							else
								t1 = floor(t1 + 0.5);
							utildavec[t] = vvec[t] = t1;
							Deltavec[t] = 0;
							if (utildavec[t] > -yvec[t])
								deltavec[t] = -1;
							else
								deltavec[t] = 1;
						}
						else {
							
							if (scsc == 1) {			
								tra1 << "@@@@@ " << jj << " " << powf(ctilda[jj], 0.5) << " " << powf(cbar, 0.5) << endl;
								scsc = 0;

								//output what is the vector found by tran_enum with the smallest index
							}
							cbar = ctilda[jj];
							for (i = jj; i <= kk; i++) {
								uvec[i] = utildavec[i];
							}
						}
					}
					else {
						t++;
						s = max(s, t);
						if (t < s) Deltavec[t] = -Deltavec[t];
						if (Deltavec[t] * deltavec[t] >= 0) Deltavec[t] += deltavec[t];
						utildavec[t] = vvec[t] + Deltavec[t];
					}
				}

				if (verb) {
					tt1 = GetTime() - tt1;
					enum_time += tt1;
				}

				NumIterations++;



			}
			else { //orth_enum

					kk = m;
					cs = ncs;
					cbar = c[jj];

					s = kk - cs;
					t = kk - 1;
					stch[0] = 3;
					stch[1] = 0;				
					qtch = 0;
					for (i = 0; i<stch[0]; i++) {
						// set initial values for x_n

						for (j = jj; j <= kk; j++) {

							if (j != kk)
								tch[qtch][i][j] = 0;
							else if(i==0){
								tch[qtch][i][j] = i;
							}
							else if(i==1){
								tch[qtch][i][j] = floor(powf(cbar / (kk - jj + 1) / c[kk], 0.5));
								if(tch[qtch][i][j]==0) // if initial values exist repetition, an extra addition is called.
									tch[qtch][i][j]=2;
							}
							else if(i==2){
								tch[qtch][i][j] = ceil(powf(cbar / (kk - jj + 1) / c[kk], 0.5));
							}
						}
						tch[qtch][i][0] = 0;
						tch[qtch][i][kk + 1] = i*i*c[kk];
					}
	

					for (i = jj; i <= kk; i++) {
						sn[i] = c[i];
					}

					for (iii = 0; iii < stch[qtch]; iii++) {
						t = s + cs - 1;
						for (j = jj; j <= kk; j++) {
							uvec[j] = tch[qtch][iii][j];
						}

						ylen[t + 1] = tch[qtch][iii][kk + 1];

						for (i = jj; i <= t; i++) {

							un[t + 1][i] = 0;
							for (j = t + 1; j <= kk; j++) {

								un[t + 1][i] += uvec[j] * mu[j][i];
							}
							
						}
						//store initial values


						for (i = s; i <= t; i++) {
							des[i] = powf(cbar / (kk - jj + 1) / c[i], 0.5);
							jw[i] = 0;
						}
						//compute average values

						for (t = kk - 1; t >=jj; t--) {//a depth first search

							if (t < s) { //when t<s, x_i has the unique choice
								if (jw[t] == 0) {
									ts1 = - un[t + 1][t];

									if (ts1 >= 0)
										ts1 = ceil(ts1 - 0.5);
									else
										ts1 = floor(ts1 + 0.5);
									xz[int(jw[t])][t] = ts1;
									jw[t]++;
								}
								deltavec[t] = 0;
							}
							else {
								
								if (jw[t] == 0) {
									//when t>=s, if x_i's chosens are not provided, compute them.

									ts1 = (des[t]) - un[t + 1][t];
									if (ts1 >= 0)
										ts1 = ceil(ts1 - 0.5);
									else
										ts1 = floor(ts1 + 0.5);

									if (des[t] != 0) {

										ts2 = -1 * (des[t]) - un[t + 1][t];
										if (ts2 >= 0)
											ts2 = ceil(ts2 - 0.5);
										else
											ts2 = floor(ts2 + 0.5);
									}
									else {
										ts2 = ts1;
									}


									ts5 = -un[t + 1][t];
									if (ts5 >= 0)
										ts5 = ceil(ts5 - 0.5);
									else
										ts5 = floor(ts5 + 0.5);

									//the three chosen are offerred


									if ((abs(abs(ts1 - ((des[t]) - un[t + 1][t])) - 0.5) < 0.1)||(abs(abs(ts2 - (-1 * (des[t]) - un[t + 1][t])) - 0.5) < 0.1)) {//not close enough, other chosen will be added					

										ts3 = ts1;
										if ((abs(abs(ts1 - ((des[t]) - un[t + 1][t])) - 0.5) < 0.1)) {

											if ((des[t]) - un[t + 1][t] > ts1) {
												ts3 = ts1;
												ts1++;
											}
											else
												ts3 = ts1 - 1;
											xz[int(jw[t])][t] = ts1;
											jw[t]++;
											xz[int(jw[t])][t] = ts3;
											jw[t]++;
										}
										else
										{
											xz[int(jw[t])][t] = ts1;
											jw[t]++;

										}

										ts4 = ts2;
										if ((abs(abs(ts2 - (-1 * (des[t]) - un[t + 1][t])) - 0.5) < 0.1)) {

											if (-1 * (des[t]) - un[t + 1][t] > ts2) {
												ts4 = ts2;
												ts2++;
											}
											else
												ts4 = ts2 - 1;

											if (ts2 == ts1) {

											}
											else if (ts2 == ts3) {
												xz[int(jw[t])][t] = ts4;
												jw[t]++;
											}
											else
											{
												xz[int(jw[t])][t] = ts2;
												jw[t]++;
												xz[int(jw[t])][t] = ts4;
												jw[t]++;
											}
										}
										else if (ts2 != ts1&&ts2 != ts3)
										{
											xz[int(jw[t])][t] = ts2;
											jw[t]++;
										}

										if (ts5 != 0 && ts5 != -0) {
											pda = 1;
											for (i = 0; i < jw[t]; i++)
												if (ts5 == xz[i][t])
												{
													pda = 0;
												}

										}
										else {
											pda = 1;
											for (i = 0; i < jw[t]; i++)
												if (ts5 == xz[i][t] || -1 * ts5 == xz[i][t])
												{
													pda = 0;
												}
										}
										if (pda == 1) {
											xz[int(jw[t])][t] = ts5;
											jw[t]++;
										}

										//Duplicate removal
									}
									else {
										if (ts1 == ts2) {
											jw[t] = 1;
											xz[0][t] = ts1;
										}
										else {
											jw[t] = 2;
											if (abs(ts1 - (des[t] - un[t + 1][t])) > abs(ts2 - (des[t] - un[t + 1][t]))) {
												ts3 = ts1;
												ts1 = ts2;
												ts2 = ts3;
											}
											xz[0][t] = ts1;
											xz[1][t] = ts2;
										}
										if (ts5 != 0 && ts5 != -0) {
											pda = 1;
											for (i = 0; i < jw[t]; i++)
												if (ts5 == xz[i][t])
												{
													pda = 0;
												}
										}
										else {
											pda = 1;
											for (i = 0; i < jw[t]; i++)
												if (ts5 == xz[i][t] || -1 * ts5 == xz[i][t])
												{
													pda = 0;
												}
										}
										if (pda == 1) {
											xz[int(jw[t])][t] = ts5;
											jw[t]++;
										}

										//Duplicate removal
									}

									deltavec[t] = 0;

								}
								else {

									if (deltavec[t] < jw[t] - 1) {

										deltavec[t]++;
									}

								}
							}
							uvec[t] = xz[deltavec[t]][t];

							//compute norm
							for (i = jj; i <= t; i++) {
								if (i != t)
									tmp = uvec[t] * mu[t][i];
								else
									tmp = uvec[t];

								un[t][i] = un[t + 1][i] + tmp;
								if (i == t)
									ylen[t] = ylen[t + 1] + (un[t][i])*(un[t][i])*c[i];
							}

							if (ylen[t] < sn[t]&& ylen[t]!=0) {// store n vecotrs

								sn[t] = ylen[t];
								for (i = t; i <= kk; i++)
									snv[t][i] = uvec[i];

							}


							if (t == jj) {

								if (!(deltavec[t] < jw[t] - 1))
								{
									jw[t] = 0;
									for (i = t + 1; i < s + cs; i++) {
										if (deltavec[i] < jw[i] - 1) {
											t = i;
											t++;
											break;
										}
										else
										{
											jw[i] = 0;
										}
									}
								}
								else
									t++;
							}
							
						}//


					}//stch

					for (i = jj; i <= kk; i++)
					{
						
						if (sn[i] < c[i]) {// check successively to decide which one to update 
							tra1 <<"!!!!  "<<i<<" "<< powf(sn[i],0.5) << " " << powf(c[i],0.5) << endl;
							
							for (j = 0; j < i; j++)
								uvec[i] = 0;
							for (j = i; j <= kk; j++) {
								uvec[j] = snv[i][j];
							}
							jj = i;
							cbar = sn[i];
							break;
						}
					}

			}

			if (jj == 1)
				shuchu = 1;

			h = min(kk + 1, m);

			if ((delta - 8 * red_fudge)*c[jj] > cbar) {

				clean = 0;

				// we treat the case that the new vector is b_s (jj < s <= kk)
				// as a special case that appears to occur most of the time.

				s = 0;
				for (i = jj + 1; i <= kk; i++) {
					if (uvec[i] != 0) {
						if (s == 0)
							s = i;
						else
							s = -1;
					}
				}
				if (s == 0) Error("BKZ_FP: internal error");

				if (s > 0) {
					// special case

					NumTrivial++;

					for (i = s; i > jj; i--) {
						// swap i, i-1
						swap(B(i - 1), B(i));
						if (U) swap((*U)(i - 1), (*U)(i));
						tp = B1[i - 1]; B1[i - 1] = B1[i]; B1[i] = tp;
						t1 = b[i - 1]; b[i - 1] = b[i]; b[i] = t1;
					}

					// cerr << "special case\n";
					new_m = ll_LLL_FP(B, U, delta, 0, check,
						B1, mu, b, c, h, jj, quit);

					if (new_m != h) Error("BKZ_FP: internal error");
					if (quit) break;
				}
				else {
					// the general case

					NumNonTrivial++;

					for (i = 1; i <= n; i++) conv(B(m + 1, i), 0);

					if (U) {
						for (i = 1; i <= m_orig; i++)
							conv((*U)(m + 1, i), 0);
					}

					for (i = jj; i <= kk; i++) {
						if (uvec[i] == 0) continue;
						conv(MU, uvec[i]);
						RowTransform2(B(m + 1), B(i), MU);
						if (U) RowTransform2((*U)(m + 1), (*U)(i), MU);
					}

					for (i = m + 1; i >= jj + 1; i--) {
						// swap i, i-1
						swap(B(i - 1), B(i));
						if (U) swap((*U)(i - 1), (*U)(i));
						tp = B1[i - 1]; B1[i - 1] = B1[i]; B1[i] = tp;
						t1 = b[i - 1]; b[i - 1] = b[i]; b[i] = t1;
					}

					for (i = 1; i <= n; i++) {
						conv(B1[jj][i], B(jj, i));
						CheckFinite(&B1[jj][i]);
					}

					b[jj] = InnerProduct(B1[jj], B1[jj], n);
					CheckFinite(&b[jj]);

					if (b[jj] == 0) Error("BKZ_FP: internal error");

					// remove linear dependencies

					// cerr << "general case\n";

					new_m = ll_LLL_FP(B, U, delta, 0, 0, B1, mu, b, c, kk + 1, jj, quit);

					if (new_m != kk) Error("BKZ_FP: internal error");

					// remove zero vector

					for (i = kk + 2; i <= m + 1; i++) {
						// swap i, i-1
						swap(B(i - 1), B(i));
						if (U) swap((*U)(i - 1), (*U)(i));
						tp = B1[i - 1]; B1[i - 1] = B1[i]; B1[i] = tp;
						t1 = b[i - 1]; b[i - 1] = b[i]; b[i] = t1;
					}

					quit = 0;
					if (check) {
						for (i = 1; i <= kk; i++)
							if ((*check)(B(i))) {
								quit = 1;
								break;
							}
					}

					if (quit) break;

					if (h > kk) {
						// extend reduced basis

						new_m = ll_LLL_FP(B, U, delta, 0, check,
							B1, mu, b, c, h, h, quit);

						if (new_m != h) Error("BKZ_FP: internal error");
						if (quit) break;
					}
				}

				z = 0;
			}
			else {
				// LLL_FP
				// cerr << "progress\n";
				NumNoOps++;

				if (!clean) {
					new_m =
						ll_LLL_FP(B, U, delta, 0, check, B1, mu, b, c, h, h, quit);
					if (new_m != h) Error("BKZ_FP: internal error");
					if (quit) break;
				}

				z++;
			}

			if (shuchu == 1) {// output current basis
		
				if (ol.is_open())
					ol.close();
				
				ol.open(shc);
				comp_sum = 0;
				for (i = 1; i <= m; i++) {
					ol << B[i - 1] << " " << powf(b[i], 0.5) << endl;
					comp_sum += powf(b[i], 0.5);
				}
				comp_st2 = GetTime();
				comp << hl << " " << beta << " " << powf(b[1], 0.5) << " " << comp_sum << " " << comp_st2 - comp_st1 << endl;
			}

		}
	}


	if (verb) {
		BKZStatus(GetTime(), enum_time, NumIterations, NumTrivial, NumNonTrivial,
			NumNoOps, m, B);
	}

	// clean up


	if (m_orig > m) {
		// for consistency, we move zero vectors to the front

		for (i = m + 1; i <= m_orig; i++) {
			swap(B(i), B(i + 1));
			if (U) swap((*U)(i), (*U)(i + 1));
		}

		for (i = 0; i < m; i++) {
			swap(B(m_orig - i), B(m - i));
			if (U) swap((*U)(m_orig - i), (*U)(m - i));
		}
	}

	B.SetDims(m_orig, n);
	BB = B;

	if (U) {
		U->SetDims(m_orig, m_orig);
		*UU = *U;
	}

	for (i = 1; i <= m_orig + 1; i++) {
		delete[] B1[i];
	}

	delete[] B1;

	for (i = 1; i <= m_orig + 1; i++) {
		delete[] mu[i];
	}

	delete[] mu;

	delete[] c;
	delete[] b;
	delete[] ctilda;
	delete[] vvec;
	delete[] yvec;
	delete[] uvec;
	delete[] utildavec;
	delete[] Deltavec;
	delete[] deltavec;

	return m;
}







long BKZ_FP_MIX(mat_ZZ& BB, double delta,
	long beta, long prune, LLLCheckFct check, long verb, int cs, char* qr2)
{
	verbose = verb;
	RR_GS_time = 0;
	NumSwaps = 0;
	if (verbose) {
		StartTime = GetTime();
		LastTime = StartTime;
	}

	if (delta < 0.50 || delta >= 1) Error("BKZ_FP: bad delta");
	if (beta < 2) Error("BKZ_FP: bad block size");

	return MBKZ(BB, 0, delta, beta, prune, check, cs, qr2);
}

NTL_END_IMPL
