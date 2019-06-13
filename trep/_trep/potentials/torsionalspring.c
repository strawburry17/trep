#include <Python.h>
#include "structmember.h"
#define TREP_MODULE
#include "../trep.h"

typedef struct
{
    Potential potential; // Inherits from Potential
    Frame *frame1;
    Frame *frame2;
    double k;
    double q0;
} TorsionalSpringPotential;



static double V(TorsionalSpringPotential *self)
{
    double theta;
    mat4x4 res;

    mul_mm4(res,
    	*Frame_g_inv(self->frame1),
    	*Frame_g(self->frame2));

    if(res[0][0]==1)
        theta = acos(res[1][1]);	
    else if(res[1][1]==1||res[2][2]==1)
        theta = acos(res[0][0]);
    else{
    	printf("frames are not coaxial");
    	return;
    }

    return 0.5 * self->k * (theta - self->q0) * (theta - self->q0);
}

static double V_dq(TorsionalSpringPotential *self, Config *q1)
{
    double theta, dtheta;
    mat4x4 mat1;
    mat4x4 mat2;
    mat4x4 mat1_dq;
    mat4x4 mat2_dq;
    mat4x4 res;

    mul_mm4(res,
	*Frame_g_inv(self->frame1),
	*Frame_g(self->frame2));

    copy_mat4x4(mat1,
    	*Frame_g_inv(self->frame1));
    copy_mat4x4(mat2,
    	*Frame_g(self->frame2));
    copy_mat4x4(mat1_dq,
    	*Frame_g_inv_dq(self->frame1, q1));
    copy_mat4x4(mat2_dq,
    	*Frame_g_dq(self->frame2, q1));
  
    if(res[0][0]==1){
        theta = acos(res[1][1]);
		dtheta= -1/sqrt(1-res[1][1]*res[1][1])*(mat1_dq[1][0]*mat2[0][1]+mat1[1][0]*mat2_dq[0][1]+mat1_dq[1][1]*mat2[1][1]+mat1[1][1]*mat2_dq[1][1]+mat1_dq[1][2]*mat2[2][1]+mat1[1][2]*mat2_dq[2][1]);
    }
    if(res[1][1]==1||res[2][2]==1){
        theta = acos(res[0][0]);
		dtheta= -1/sqrt(1-res[0][0]*res[0][0])*(mat1_dq[0][0]*mat2[0][0]+mat1[0][0]*mat2_dq[0][0]+mat1_dq[0][1]*mat2[1][0]+mat1[0][1]*mat2_dq[1][0]+mat1_dq[0][2]*mat2[2][0]+mat1[0][2]*mat2_dq[2][0]);
    }

    printf("%lf\t%lf\t%lf\n",theta,dtheta,self->q0);
    
    if (isnan(dtheta) && (self->q0 == 0))
	return 0;
    return self->k * (theta - self->q0) * dtheta;

}

static double V_dqdq(TorsionalSpringPotential *self, Config *q1, Config *q2)
{
    double theta, ditheta, djtheta, ddtheta;
    mat4x4 mat1;
    mat4x4 mat2;
    mat4x4 mat1_dqi;
    mat4x4 mat2_dqi;
    mat4x4 mat1_dqj;
    mat4x4 mat2_dqj;
    mat4x4 mat1_ddq;
    mat4x4 mat2_ddq;
    mat4x4 res;

    mul_mm4(res,
	*Frame_g_inv(self->frame1),
	*Frame_g(self->frame2));

    copy_mat4x4(mat1,
    	*Frame_g_inv(self->frame1));
    copy_mat4x4(mat2,
    	*Frame_g(self->frame2));
    copy_mat4x4(mat1_dqi,
    	*Frame_g_inv_dq(self->frame1, q1));
    copy_mat4x4(mat2_dqi,
    	*Frame_g_dq(self->frame2, q1));
    copy_mat4x4(mat1_dqj,
    	*Frame_g_inv_dq(self->frame1, q2));
    copy_mat4x4(mat2_dqj,
    	*Frame_g_dq(self->frame2, q2));
    copy_mat4x4(mat1_ddq,
    	*Frame_g_inv_dqdq(self->frame1, q1,q2));
    copy_mat4x4(mat2_ddq,
    	*Frame_g_dqdq(self->frame2, q1,q2));
    if(res[0][0]==1){
        theta = acos(res[1][1]);
    	ditheta = -1/sqrt(1-res[1][1]*res[1][1])*(mat1_dqi[1][0]*mat2[0][1]+mat1[1][0]*mat2_dqi[0][1]+mat1_dqi[1][1]*mat2[1][1]+mat1[1][1]*mat2_dqi[1][1]+mat1_dqi[1][2]*mat2[2][1]+mat1[1][2]*mat2_dqi[2][1]);
    	djtheta = -1/sqrt(1-res[1][1]*res[1][1])*(mat1_dqj[1][0]*mat2[0][1]+mat1[1][0]*mat2_dqj[0][1]+mat1_dqj[1][1]*mat2[1][1]+mat1[1][1]*mat2_dqj[1][1]+mat1_dqj[1][2]*mat2[2][1]+mat1[1][2]*mat2_dqj[2][1]);
    	ddtheta = -res[1][1]*(mat1_dqj[1][0]*mat2[0][1]+mat1[1][0]*mat2_dqj[0][1]+mat1_dqj[1][1]*mat2[1][1]+mat1[1][1]*mat2_dqj[1][1]+mat1_dqj[1][2]*mat2[2][1]+mat1[1][2]*mat2_dqj[2][1])/pow((1-res[1][1]*res[1][1]),1.5)*
    	(mat1_dqi[1][0]*mat2[0][1]+mat1[1][0]*mat2_dqi[0][1]+mat1_dqi[1][1]*mat2[1][1]+mat1[1][1]*mat2_dqi[1][1]+mat1_dqi[1][2]*mat2[2][1]+mat1[1][2]*mat2_dqi[2][1])-
    	1/sqrt(1-res[1][1]*res[1][1])*(mat1_ddq[1][0]*mat2[0][1]+mat1_dqj[1][0]*mat2_dqi[0][1]+mat1_ddq[1][1]*mat2[1][1]+mat1_dqj[1][1]*mat2_dqi[1][1]+mat1_ddq[1][2]*mat2[2][1]+mat1_dqj[1][2]*mat2_dqi[2][1])-
    	1/sqrt(1-res[1][1]*res[1][1])*(mat1_dqi[1][0]*mat2_dqj[0][1]+mat1[1][0]*mat2_ddq[0][1]+mat1_dqi[1][1]*mat2_dqj[1][1]+mat1[1][1]*mat2_ddq[1][1]+mat1_dqi[1][2]*mat2_dqj[2][1]+mat1[1][2]*mat2_ddq[2][1]);
    }
    if(res[1][1]==1||res[2][2]==1){
        theta = acos(res[0][0]);
    	ditheta = -1/sqrt(1-res[0][0]*res[0][0])*(mat1_dqi[1][0]*mat2[0][1]+mat1[1][0]*mat2_dqi[0][1]+mat1_dqi[1][1]*mat2[1][1]+mat1[1][1]*mat2_dqi[1][1]+mat1_dqi[1][2]*mat2[2][1]+mat1[1][2]*mat2_dqi[2][1]);
    	djtheta = -1/sqrt(1-res[0][0]*res[0][0])*(mat1_dqj[1][0]*mat2[0][1]+mat1[1][0]*mat2_dqj[0][1]+mat1_dqj[1][1]*mat2[1][1]+mat1[1][1]*mat2_dqj[1][1]+mat1_dqj[1][2]*mat2[2][1]+mat1[1][2]*mat2_dqj[2][1]);
    	ddtheta = res[0][0]*(mat1_dqj[1][0]*mat2[0][1]+mat1[1][0]*mat2_dqj[0][1]+mat1_dqj[1][1]*mat2[1][1]+mat1[1][1]*mat2_dqj[1][1]+mat1_dqj[1][2]*mat2[2][1]+mat1[1][2]*mat2_dqj[2][1])/pow((1-res[0][0]*res[0][0]),1.5)*
    	(mat1_dqi[1][0]*mat2[0][1]+mat1[1][0]*mat2_dqi[0][1]+mat1_dqi[1][1]*mat2[1][1]+mat1[1][1]*mat2_dqi[1][1]+mat1_dqi[1][2]*mat2[2][1]+mat1[1][2]*mat2_dqi[2][1])-
    	1/sqrt(1-res[0][0]*res[0][0])*(mat1_ddq[1][0]*mat2[0][1]+mat1_dqj[1][0]*mat2_dqi[0][1]+mat1_ddq[1][1]*mat2[1][1]+mat1_dqj[1][1]*mat2_dqi[1][1]+mat1_ddq[1][2]*mat2[2][1]+mat1_dqj[1][2]*mat2_dqi[2][1])-
    	1/sqrt(1-res[0][0]*res[0][0])*(mat1_dqi[1][0]*mat2_dqj[0][1]+mat1[1][0]*mat2_ddq[0][1]+mat1_dqi[1][1]*mat2_dqj[1][1]+mat1[1][1]*mat2_ddq[1][1]+mat1_dqi[1][2]*mat2_dqj[2][1]+mat1[1][2]*mat2_ddq[2][1]);
    }

    return self->k * ditheta * djtheta + self->k * (theta - self->q0) * ddtheta;
 
}


static void dealloc(TorsionalSpringPotential *self)
{
    Py_CLEAR(self->frame1);
    Py_CLEAR(self->frame2);
    ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
}

static int init(TorsionalSpringPotential *self, PyObject *args, PyObject *kwds)
{
    // Note that we do not called Potential.__init__ here.  It will
    // be called by Potential.__init__.
    self->k = 0.0;
    self->q0 = 0.0;

    /* Don't touch these though */
    self->potential.V = (PotentialFunc_V)&V;
    self->potential.V_dq = (PotentialFunc_V_dq)&V_dq;
    self->potential.V_dqdq = (PotentialFunc_V_dqdq)&V_dqdq;

    return 0;
}

/* Make custom parameters available to python here. */
static PyMemberDef members_list[] = {
    {"_frame1", T_OBJECT_EX, offsetof(TorsionalSpringPotential, frame1), 0, trep_internal_doc},
    {"_frame2", T_OBJECT_EX, offsetof(TorsionalSpringPotential, frame2), 0, trep_internal_doc},
    {"_q0", T_DOUBLE, offsetof(TorsionalSpringPotential, q0), 0, trep_internal_doc},
    {"_k", T_DOUBLE, offsetof(TorsionalSpringPotential, k), 0, trep_internal_doc},
    {NULL}  /* Sentinel */
};


extern PyTypeObject PotentialType;
PyTypeObject TorsionalSpringPotentialType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_trep._TorsionalSpringPotential", /*tp_name*/
    sizeof(TorsionalSpringPotential),  /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)dealloc,       /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    0,                         /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,                         /* tp_methods */
    members_list,              /* tp_members */
    0,                         /* tp_getset */
    &PotentialType,            /* tp_base */   
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)init,            /* tp_init */
    0,                         /* tp_alloc */
    0,                       /* tp_new */
};


