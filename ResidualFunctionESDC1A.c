/*
 * Module: ESDC1A Residual Function 
 * --------------------
 * 
 * Description: builds the residual function equations for the ESDC1A model.
 *
 */
//
#include "../../include/TS3ph.h"
#undef __FUNCT__
#define __FUNCT__ "ResidualFunctionESDC1A"
PetscErrorCode ResidualFunctionESDC1A(TS3phDynCtx* Gen, TS3phNetCtx* Net, PetscInt i)
{

  PetscErrorCode ierr;
  PetscScalar    *xgen_arr, *fgen_arr, *v_arr, *ugen_arr;
  PetscScalar    *xdot_arr;
  PetscScalar    *vgen, theta_a, theta_b, theta_c;
  PetscScalar    Vc, Vll, Vr, Efd, Rf;
  PetscScalar    Vc_dot, Vll_dot, Vr_dot, Efd_dot, Rf_dot;
  PetscInt       Vc_idx, Vll_idx, Vr_idx, Efd_idx, Rf_idx;
  PetscScalar    Tr, Ka, Ta, Tb, Tc, Vrmax, Vrmin;
  PetscScalar    Ke, Te, Kf, Tf1, E1, SeE1, E2, SeE2;
  PetscScalar    Vref, Ec, Vothsg, Voel, Vuel;
  PetscScalar    fgen_Vr, fgen_Efd;

  PetscInt       state_index, idx_exc;
  Vec            Xgen = Gen->xdyn, Fgen = Gen->f_dyn, V = Net->V;
  Vec            Xdot = Gen->xdotdyn, Ugen = Gen->udyn;

  GENDATA        *Gdata = Gen->GenData;
  ESDC1A_MODEL   *ESDC1A = Gdata->EXC_MODEL->ESDC1A;
  MAC_INFO_data  *MAC_INFO = Gdata->MAC_INFO;
  PetscInt       bus_index, seq_number;

  /*External vectors*/
  PetscFunctionBegin;
  ierr = VecGetArray(Xgen,&xgen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Xdot,&xdot_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Ugen,&ugen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Fgen,&fgen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(V,&v_arr);CHKERRQ(ierr);

  idx_exc = MAC_INFO->exc_idx[i] ; //TODO: This has to be addressed ASAP

  /* State variable indexes */

  Vc_idx      = ESDC1A->VC_idx[idx_exc];
  Vll_idx     = ESDC1A->VLL_idx[idx_exc];
  Vr_idx      = ESDC1A->VR_idx[idx_exc];
  Efd_idx     = ESDC1A->EFD_idx[idx_exc];
  Rf_idx      = ESDC1A->VF_idx[idx_exc];

  /*State Variables*/
  state_index = ESDC1A->VC_idx[idx_exc];
  Vc          = xgen_arr[Vc_idx];
  Vll         = xgen_arr[Vll_idx];
  Vr          = xgen_arr[Vr_idx];
  Efd         = xgen_arr[Efd_idx];
  Rf          = xgen_arr[Rf_idx];

  /*Derivatives state variables*/
  Vc_dot      = xdot_arr[Vc_idx];
  Vll_dot     = xdot_arr[Vll_idx];
  Vr_dot      = xdot_arr[Vr_idx];
  Efd_dot     = xdot_arr[Efd_idx];
  Rf_dot      = xdot_arr[Rf_idx];

  /* Parameters*/
  Tr    = ESDC1A->Tr[idx_exc];
  Ka    = ESDC1A->Ka[idx_exc];
  Ta    = ESDC1A->Ta[idx_exc];
  Tb    = ESDC1A->Tb[idx_exc];
  Tc    = ESDC1A->Tc[idx_exc];
  Vrmax = ESDC1A->Vrmax[idx_exc];
  Vrmin = ESDC1A->Vrmin[idx_exc];
  Ke    = ESDC1A->Ke[idx_exc];
  Te    = ESDC1A->Te[idx_exc];
  Kf    = ESDC1A->Kf[idx_exc];
  Tf1   = ESDC1A->Tf1[idx_exc];
  E1    = ESDC1A->E1[idx_exc];
  SeE1  = ESDC1A->SeE1[idx_exc];
  E2    = ESDC1A->E2[idx_exc];
  SeE2  = ESDC1A->SeE2[idx_exc];

  /*External and control variables*/

  Vref = ugen_arr[2*i];

  if (ESDC1A->VOTHSG_idx[idx_exc] != -1){  /* Stabilizer */
    Vothsg = xgen_arr[ESDC1A->VOTHSG_idx[idx_exc]];
  } else {
    Vothsg = 0;
  }


  if (ESDC1A->VUEL_idx[idx_exc] != -1){  /* Stabilizer */
    Vuel   = xgen_arr[ESDC1A->VUEL_idx[idx_exc]];
  } else {
    Vuel   = -100;
  }


  if (ESDC1A->VOEL_idx[idx_exc] != -1){  /* Stabilizer */
    Voel   = xgen_arr[ESDC1A->VOEL_idx[idx_exc]];
  } else {
    Voel   = 0;
  }


  bus_index  = MAC_INFO->gen_bus[i];
  seq_number = bus_index;

  if (ESDC1A->ECOMP_idx[idx_exc] != -1){  /* Compensator */
    Ec = xgen_arr[ESDC1A->ECOMP_idx[idx_exc]];
  } else { /* Input is bus voltage*/
    vgen = v_arr + twonphase*(seq_number);
    Ec = sqrt(vgen[0]*vgen[0] + vgen[1]*vgen[1] + 
    vgen[2]*vgen[2] + vgen[3]*vgen[3] + vgen[4]*vgen[4] + 
    vgen[5]*vgen[5])/sqrt(3);
  }

  /*Function mismatch construction*/

  PetscScalar Vf;
  Vf = (Kf/Tf1)*Efd - Rf;

  /*Vc*/

  if(Tr > 0){
    switch (Gen->net_mode) {

      case 0:
        fgen_arr[Vc_idx] = (1/Tr)*(Ec-Vc) - Vc_dot;
        break;
      case 1:
        fgen_arr[Vc_idx] = (1/Tr)*(Ec-Vc);
        break;
      case 2:
        fgen_arr[Vc_idx] = 0;
        break;
      default:
        printf("Error ResidualFunctionESDC1A\n");
        break;

    }

  }else{ // Algebraic                                                                                                              
    fgen_arr[Vc_idx] = Ec - Vc;
  }

  /*Vll*/

  PetscScalar aux_vll;
  aux_vll = Vref - Vc + Vothsg + Voel - Vf;

  if(Tb > 0){

    switch (Gen->net_mode) {

      case 0:
        fgen_arr[Vll_idx] = (1/Tb)*((1-(Tc/Tb))*aux_vll - Vll) - Vll_dot;
        break;
      case 1:
        fgen_arr[Vll_idx] = (1/Tb)*((1-(Tc/Tb))*aux_vll - Vll);
        break;
      case 2:
        fgen_arr[Vll_idx] = 0;
        break;
      default:
        printf("Error ResidualFunctionESDC1A\n");
        break;

    }

  }else{
    fgen_arr[Vll_idx] = aux_vll - Vll;
  }

  /* Vr */

  PetscScalar HV, Vr_aux;

  if (Tb == 0) {
    Vr_aux = Vll;
  }else{
    Vr_aux = (Tc/Tb)*aux_vll + Vll;
  }

  if (Vr_aux > Vuel) { // High Value Gate                                                                                          
    HV = Vr_aux;
  }else{
    HV = Vuel;
  }

  if(Ta > 0){

    switch (Gen->net_mode) {

      case 0:
        fgen_Vr = (1/Ta)*(Ka*HV - Vr);
        fgen_arr[Vr_idx] = fgen_Vr - Vr_dot;
        break;
      case 1:
        fgen_Vr = (1/Ta)*(Ka*HV - Vr);
        fgen_arr[Vr_idx] = fgen_Vr;
        break;
      case 2:
        fgen_Vr = 0;
        fgen_arr[Vr_idx] = fgen_Vr;
        break;
      default:
        printf("Error ResidualFunctionESDC1A\n");
        break;

    }

    /* Vr limit enforcement */
    ESDC1A->VR_lflag[idx_exc] = 0; /*Flag limit deactivated by default*/
    if((Vr >= Vrmax - 2.0e-16) && (fgen_Vr > 0)){
      fgen_arr[Vr_idx] = Vr - Vrmax;
      ESDC1A->VR_lflag[idx_exc] = 1;
    }else if ((Vr <= Vrmin + 2.0e-16) && (fgen_Vr < 0)){
      fgen_arr[Vr_idx] = Vr - Vrmin;
      ESDC1A->VR_lflag[idx_exc] = -1;
    }

  }else{
    fgen_arr[Vr_idx] = Ka*HV - Vr;

    /* Vr limit enforcement */
    ESDC1A->VR_lflag[idx_exc] = 0; /*Flag limit deactivated by default*/
    if(Vr >= Vrmax - 2.0e-16){
      fgen_arr[Vr_idx] = Vr - Vrmax;
      ESDC1A->VR_lflag[idx_exc] = 1;
    }else if (Vr <= Vrmin + 2.0e-16){
      fgen_arr[Vr_idx] = Vr - Vrmin;
      ESDC1A->VR_lflag[idx_exc] = -1;
    }

  }


  /*Efd*/

  /*first we calculate saturation parameters*/
  /*Suggestion: This should be done in the initialization*/
  PetscScalar SeEFD, A, B, Z, Z1, A1, B1;


  if (Gen->net_mode != 2) {

    if ((SeE2 != 0) && (E2 != 0)){ /*check for saturation data correctness*/
      Z = sqrt((SeE1*E1)/(SeE2*E2));
      A = (Z*E2-E1)/(Z-1);
      B = pow((sqrt(SeE1*E1)-sqrt(SeE2*E2))/(E1-E2), 2);
      if (Efd >= A){
        fgen_Efd = (1/Te)*(Vr - Ke*Efd - B*pow(Efd-A,2));
      }else{
        fgen_Efd = (1/Te)*(Vr - Ke*Efd);
      }
      fgen_arr[Efd_idx] = fgen_Efd;    
    }else {

      fgen_Efd  = (1/Te)*(Vr - Ke*Efd);
      fgen_arr[Efd_idx] = fgen_Efd;        
    }

    if (Gen->net_mode == 0)  fgen_arr[Efd_idx] -= Efd_dot;

    ESDC1A->EFD_lflag[idx_exc] = 0; /*Flag limit deactivated by default*/

    if((Efd <= 0) && (fgen_Efd < 0)){
      fgen_arr[Efd_idx] = Efd;
      ESDC1A->EFD_lflag[idx_exc] = 1; /*Activate flag limit*/
    }

  }else{
    fgen_arr[Efd_idx] = 0;
  }




  /*Rf*/

  switch (Gen->net_mode) {

    case 0:
      fgen_arr[Rf_idx] = (1/Tf1)*Vf - Rf_dot;
      break;
    case 1:
      fgen_arr[Rf_idx] = (1/Tf1)*Vf;
      break;
    case 2:
      fgen_arr[Rf_idx] = 0;
      break;
    default:
      printf("Error ResidualFunctionESDC1A\n");
      break;

  }

  /*Restore arrays*/

  ierr = VecRestoreArray(Xdot,&xdot_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xgen,&xgen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Fgen,&fgen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ugen,&ugen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(V,&v_arr);CHKERRQ(ierr);
}
