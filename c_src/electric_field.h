/* This file is part of SOMA.

 SOMA is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SOMA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SOMA_ELECTRIC_FIELD_H
#define SOMA_ELECTRIC_FIELD_H

#include "soma_config.h"
#include "phase.h"


//! Top level struct for electric field implementation.
typedef struct ElectricField{
	soma_scalar_t *eps;				    //!< Array containing the dielectric constants for all particle types.
    soma_scalar_t *eps_arr;             //!< Array that saves information of calculated dielectric constant field.
    uint8_t *electrodes;                //!< Array containing the electrode positions.
    uint64_t iter_limit;                //!< Value that determines the upper limit of iterations to solve the electric field.
    soma_scalar_t thresh_iter;          //!< Value that determines the threshold to stop iterative solution of the electric field.
    uint64_t amt_iter;                  //!< Value that stores the amount of iterations to solve the electric field (per MC step).
    soma_scalar_t *Epot;                //!< Array containing the electric potential field.
    soma_scalar_t *Epot_tmp;            //!< Temporary array that contains the electric potential field after MC step.
    soma_scalar_t *pre_deriv;           //!< Array containing precomputed derivatives of the dielectric constant field.
    soma_scalar_t *H_el_field;          //!< Array containing the cell-wise contribution to electrostatic energy hamiltonian.
    soma_scalar_t *E_field;             //!< Array containing the electric field.
    soma_scalar_t H_el;                 //!< Electrostatic energy hamiltonia.
    soma_scalar_t *omega_field_el;      //!< Array containing electrotatic energy contribution to omega fields.
    soma_scalar_t  sqrt_Nbar;           //!< Value of \sqrt{\hat{N}}.
    uint8_t stride;                     //!< Value of stride used for convolution.
    soma_scalar_t *kernel;              //!< Array containing kernel for convolution.
    soma_scalar_t *kernel_norm_field;   //!< Array containing kernel normalization field.
    soma_scalar_t *kernel_blur;         //!< Array containing simple blur kernel to smooth Epot after deconvolution.
    uint8_t kernel_dim;                 //!< Value containing kernel dimension.
    uint8_t kernel_rad;                 //!< Value of kernel radius (integer).
    soma_scalar_t kernel_sigma;         //!< Value for standard deviation (sigam) of gaussian kernel.
    unsigned int conv_nx;               //!< Convoluted x-spatial discretization.
    unsigned int conv_ny;               //!< Convoluted y-spatial discretization.
    unsigned int conv_nz;               //!< Convoluted z-spatial discretization.
    uint64_t n_cells_conv;              //!< Number of cells after convolution.
    bool el_pos_xy;                     //!< Bool to determine if electrode position is in xy-plane.
    bool el_pos_xz;                     //!< Bool to determine if electrode position is in xz-plane.
    bool el_pos_yz;                     //!< Bool to determine if electrode position is in yz-plane.
    uint8_t x_offset;                   //!< Offset for convoluted x axis; used to skip planes exhibiting electrodes during iteration.
    uint8_t y_offset;                   //!< Offset for convoluted y axis; used to skip planes exhibiting electrodes during iteration.
    uint8_t z_offset;                   //!< Offset for convoluted z axis; used to skip planes exhibiting electrodes during iteration.
    soma_scalar_t *electrodes_conv;     //!< Array containing the convoluted values of the electrode field.
    soma_scalar_t *eps_arr_conv;        //!< Array containing the convoluted values of the dielectric constant field.
    soma_scalar_t *pre_deriv_conv;      //!< Array containing the precomputed derivatives of the convoluted dielectric constant field.
    soma_scalar_t *Epot_conv;           //!< Array containing the convoluted electric potential field.
    soma_scalar_t *Epot_tmp_conv;       //!< Temporary array containing the convoluted electric potential field after MC step.
    
} ElectricField;

/*! Helper function to read the electric field array, electrode array and dielectric constants from the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \returns Errorcode */
int read_electric_field_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to write the electric field array, electrode array and dielectric constants to the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \returns Errorcode */
int write_electric_field_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to calculate sqrt(N_bar)
    \private
    \param p Fully CPU initialized Phase struct
    \returns Errorcode */
void calc_sqrt_Nbar(struct Phase *const p);

/*! Helper function to determine electrode locations in order to resolve non-periodic boundaries
    \private
    \param p Fully CPU initialized Phase struct
    \returns Errorcode */
int init_efield(struct Phase *const p);

/*! Helper function to allocate and compute gaussian kernel for convolution
    \private
    \param p Fully CPU initialized Phase struct
    \returns Errorcode */
int init_kernel(struct Phase *const p);

/*! Helper function to determine and allocate arrays () for convolution
    \private
    \param p Fully CPU initialized Phase struct
    \returns Errorcode */
int init_convolution(struct Phase *const p);

/*! Helper function to copy the ef data to the device
    \private
    \param p Fully CPU initialized Phase struct
    \returns Errorcode */
int copyin_electric_field(struct Phase *p);

/*! Helper function delete the ef data from the device and copy it to the CPU memory
    \private
    \param p Fully CPU initialized Phase struct
    \returns Errorcode */
int copyout_electric_field(struct Phase *p);

/*! Helper function to update the host with the ef data
    \private
    \param p Fully initialized Phase struct
    \returns Errorcode */
int update_self_electric_field(const struct Phase *const p);

// /*! Helper function to compute the value of sqrt{\bar{N}} defined as the amount of polymer chains per volume R_e^3
//    \private
//    \param p Phase describing the system */
// void calc_sqrt_Nbar(struct Phase *const p);

/*! Helper function to compute the dielectric field from the densities of individual particle types (welling2014, eq. 83)
    \private
    \param p Phase describing the system
    \returns sqrt{\bar{N}} as uint */
void calc_dielectric_field(struct Phase *const p);

/*! Helper funtion to precompute derivatives of the dielectric field (welling2014, eq. 85)
    \private
    \param p Phase describing the system */
void pre_derivatives(struct Phase *const p);

/*! Helper function to iterate solution for derivatives of the electric potential field 
 *  (using precomputed derivatives of the dielectric field)
    \private
    \param p Phase describing the system
    \returns maximum of dE_pot/dr */
soma_scalar_t iterate_field(struct Phase *const p);

/*! Helper function to convolute the dielectric field array using the gauss kernel
    \private
    \param p Phase describing the system */
void convolution_eps_arr(struct Phase *const p);

/*! Helper funtion to precompute derivatives of the convoluted dielectric field (welling2014, eq. 85)
    \private
    \param p Phase describing the system */
void pre_derivatives_conv(struct Phase *const p);

/*! Helper function to iterate solution for derivatives of the convoluted electric potential field 
 *  (using precomputed derivatives of the convoluteddielectric field)
    \private
    \param p Phase describing the system
    \returns maximum of dE_pot/dr */
soma_scalar_t iterate_field_conv(struct Phase *const p);

/*! Helper function to deconvolute the electric potential field array to its original size
    \private
    \param p Phase describing the system */
void deconvolution_Epot(struct Phase *const p);

/*! Main routine, calculates electrostatic energy contribution per cell and total (welling2017, eq. 7)
    \private
    \param p Phase describing the system
    \returns Errorcode */
int calc_electric_field_contr(struct Phase *const p);

/*! Helper function to free the CPU memory resources of the pc struct. The function gets automatically called by free_phase().
  \private
  \param p Initialized Phase that is in the process of deallocating its resources.
  \returns Errorcode */
int free_electric_field(struct Phase *const p);

/*! Helper function to track various parameter for debugging purposes.
  \private
  \param p Initialized Phase that is in the process of deallocating its resources.
  \param k Amount of iterations to solve electrical potential field. */
void tests(struct Phase *const p, uint64_t k);

#endif                          //SOMA_ELECTRIC_FIELD_H
