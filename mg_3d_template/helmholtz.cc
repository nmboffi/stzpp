#include <cstdlib>

#include "helmholtz.hh"

helmholtz::helmholtz(const double ax_,const double bx_,
			const double ay_,const double by_,const double az_,const double bz_,
			double k,geometry &gm,comm_buffer &com) :
	x_prd(gm.x_prd), y_prd(gm.y_prd), z_prd(gm.z_prd),
	rank(gm.rank), m(gm.m), n(gm.n), o(gm.o), mno(m*n*o),
	ai((m*gm.ip)/gm.mp), aj((n*gm.jp)/gm.np), ak((o*gm.kp)/gm.op),
	bi((m*(gm.ip+1))/gm.mp), bj((n*(gm.jp+1))/gm.np), bk((o*(gm.kp+1))/gm.op),
	sm(bi-ai), sn(bj-aj), so(bk-ak), smno(sm*sn*so),
	ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_),
	dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)), dz((bz-az)/(o-1)),
	fm(cpx(-2/(dx*dx)-2/(dy*dy)-2/(dz*dz),k*k)), fm_inv(1./fm),
	fm_full(-2/(dx*dx)-2/(dy*dy)-2/(dz*dz)+k*k), fex(1./(dx*dx)),
	fey(1./(dy*dy)), fez(1./(dz*dz)), b(new cpx[11*smno]), x(b+smno),
	mg(*this,gm,com) {

	// Set the stencil entries
	const cpx stencil_[27]={0,0,0,0,fez,0,0,0,0,
			                0,fey,0,fex,fm,fex,0,fey,0,
							0,0,0,0,fez,0,0,0,0};
	memcpy(stencil,stencil_,27*sizeof(cpx));
	*fix=fm;

	// Set required constants from the top multigrid region
	region<cpx,cpx>* &rp=mg.reg[0];
	reg_r0=rp->r0;
	reg_x0=rp->x0;
	rmg=rp->mg;
	rng=rp->ng;
}

/** The class destructor frees the dynamically allocated memory. */
helmholtz::~helmholtz() {
	delete [] b;
}

/** Solves the Helmholtz problem using the preconditioned Bi-CGSTAB
 * method. */
void helmholtz::solve() {
	int ij,l=0;
	cpx rho=1,alpha=1,omega=1,nrho,beta;
	cpx *p=b+2*smno,*v=b+3*smno,*rhat=b+4*smno,*r=b+5*smno,*s=b+6*smno,
	    *zz=b+7*smno,*kt=b+8*smno,*y=b+9*smno,*t=b+10*smno;
	double l2;

	// Set up the multigrid hierarchy
	mg.setup_matrices(*this);

	// Set up
	for(ij=0;ij<smno;ij++) v[ij]=p[ij]=0;

	// Choose the arbitrary vector rhat to just be equal to the initial r
	mul_fa(x,r);
	for(ij=0;ij<smno;ij++) rhat[ij]=r[ij]=b[ij]-r[ij];

	// Do the Bi-CGSTAB iteration
	l2=l2_error(r);
	if(rank==0) printf("Iter %d, L2 error %.10g\n",l,l2);
	while(true) {
		nrho=iprod(rhat,r);
		beta=(nrho/rho)*(alpha/omega);
		rho=nrho;
		for(ij=0;ij<smno;ij++) p[ij]=r[ij]+beta*(p[ij]-omega*v[ij]);
		pc_solve(p,y);
		mul_fa(y,v);
		alpha=rho/iprod(rhat,v);
		for(ij=0;ij<smno;ij++) s[ij]=r[ij]-alpha*v[ij];
		pc_solve(s,zz);
		mul_fa(zz,t);
		pc_solve(t,kt);
		omega=iprod(kt,zz)/iprod(kt,kt);
		for(ij=0;ij<smno;ij++) x[ij]=x[ij]+alpha*y[ij]+omega*zz[ij];
		l++;
		for(ij=0;ij<smno;ij++) r[ij]=s[ij]-omega*t[ij];
		l2=l2_error(r);
		if(rank==0) printf("Iter %d, L2 error %.10g\n",l,l2);
		if(l==20||l2<1e-20) return;
	}
}

/** Calculates the inner product of two complex fields.
 * \param[in] (u,v) pointers to the two fields. */
cpx helmholtz::iprod(cpx *u,cpx *v) {
	cpx w=std::conj(*u)*(*v);
	for(int ij=1;ij<smno;ij++) w+=std::conj(u[ij])*v[ij];

	// Sum the contributions from all processors
	double q[4];
	*q=w.real();q[1]=w.imag();
	MPI_Allreduce(q,q+2,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	return cpx(q[2],q[3]);
}

/** Solves the preconditioning problem, to test the multigrid code. */
void helmholtz::pc_test() {
	int i,j,k;
	double l2;
	mg.setup_matrices(*this);

	memcpy(reg_r0,b,smno*sizeof(cpx));

	for(i=0;i<10;i++) {
		l2=mg.l2_error();
		if(rank==0) printf("MG iter %d, L2 error %.10g\n",i,l2);
		mg.v_cycle();
	}
	l2=mg.l2_error();
	if(rank==0) printf("MG iter %d, L2 error %.10g\n",i,l2);

	cpx *xp=x;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++)
		*(xp++)=reg_x0[i+rmg*(j+rng*k)];
}

/** Solves the preconditioning problem using multigrid.
 * \param[in] u a pointer to the source data.
 * \param[in] v a pointer to the solution. */
void helmholtz::pc_solve(cpx *u,cpx *v) {
	int i,j,k;
//memcpy(v,u,smno*sizeof(cpx));
//	return;

	memcpy(reg_r0,u,smno*sizeof(cpx));

	for(i=0;i<1;i++) mg.v_cycle();

	cpx *vp=v;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++)
		*(vp++)=reg_x0[i+rmg*(j+rng*k)];
}

/** Multiplies an input field by the full problem A.
 * \param[in] u a pointer to the input field to consider.
 * \param[in] v a pointer to the solution. */
void helmholtz::mul_fa(cpx *u,cpx *v) {
	int i,j,k;
	cpx *pp=u,*xp;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++)
		reg_x0[i+rmg*(j+rng*k)]=*(pp++);

	mg.reg[0]->communicate();
	pp=v;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		xp=reg_x0+i+rmg*(j+rng*k);
		*(pp++)=interior(ai+i, aj+j, ak+k)?
				*xp*fm_full+(xp[1]+xp[-1])*fex+(xp[-rmg]+xp[rmg])*fey
			    +(xp[-rmg*rng]+xp[rmg*rng])*fez:*xp*(*fix);
	}
}

/** Calculates the mean squared error of a complex field.
 * \param[in] u a pointer to the field.
 * \return The mean squared error. */
double helmholtz::l2_error(cpx *u) {
	double l2=std::norm(*u),l2s;
	for(int ij=1;ij<smno;ij++) l2+=std::norm(u[ij]);

	MPI_Allreduce(&l2,&l2s,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	return l2s/mno;
}

void helmholtz::setup_test_problem() {
	int i,j,k;
	double xx,yy,zz;

	cpx *xp=x,*bp=b;
	for(k=0;k<so;k++) {
		zz=az+dz*(k+ak);
		for(j=0;j<sn;j++) {
			yy=ay+dy*(j+aj);
			for(i=0;i<sm;i++) {
				xx=ax+dx*(i+ai);
				*(xp++)=0;
				*(bp++)=xx>0&&xx<0.3&&yy>0&&yy<0.3&&zz>0&&zz<0.3?(*fix):0;
			}
		}
	}
}

void helmholtz::output_x(const char* odir) {
	int i,j,k;
	cpx *pp=x;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++)
		reg_x0[i+rmg*(j+rng*k)]=*(pp++);

	// Output the header
	char buf[64];
	if(rank==0) {
		sprintf(buf,"%s/header",odir);
		FILE *fp=safe_fopen(buf,"w");
		fprintf(fp,"%g %g %d\n",az,bz,o);
		fclose(fp);
	}

	// Output the cross-sections in z
	for (int i = 0; i < o; i++) {
		sprintf(buf, "%s/x.%d", odir, i);
		mg.output_x(ax,bx,ay,by,buf, i, 0);
	}
}

void helmholtz::output_b(const char* odir) {
	memcpy(reg_r0,b,smno*sizeof(cpx));

	// Output the cross-sections in z
	char buf[64];
	for (int i = 0; i < o; i++) {
		sprintf(buf, "%s/r.%d", odir, i);
		mg.output_r(ax,bx,ay,by,buf, i, 0);
	}
}

/** Fills the entries in the table with stencil values for a particular gridpoint.
 * \param[in] (i,j,k) the index of the gridpoint to consider.
 * \param[in,out] en a reference to the pointer to the stencil table. */
void helmholtz::fill_entries(int i, int j, int k, cpx *&en) {
	p_fatal_error("No need to fill entries with this problem type", 1);
}

// Explicit instantiation
inline cpx mg3d_inverse(cpx a) {return 1./a;}
inline double mg3d_mod_sq(cpx a) {return std::norm(a);}
inline double mg3d_float(cpx a) {return static_cast<float>(a.real());}
#include "region.cc"
#include "multigrid.cc"
template class region<cpx,cpx>;
template class multigrid<cpx,cpx>;
