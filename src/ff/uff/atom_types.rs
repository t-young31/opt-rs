use crate::ff::uff::atom_typing::UFFAtomType;


pub(crate) const ATOM_TYPES: [UFFAtomType; 127] = [
UFFAtomType{name: "H_", atomic_symbol: "H", bridging: false, aromatic: false, valency: 1, oxidation_state: 0, r: 0.354, theta: 3.14159, x: 2.886, d: 0.105, zeta: 12.73, z_eff: 1.912}, 
UFFAtomType{name: "H_b", atomic_symbol: "H", bridging: true, aromatic: false, valency: 2, oxidation_state: 0, r: 0.460, theta: 1.45735, x: 2.886, d: 0.111, zeta: 12.0, z_eff: 1.787}, 
UFFAtomType{name: "He4+4", atomic_symbol: "He", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 0.849, theta: 3.14159, x: 3.851, d: 0.505, zeta: 11.278, z_eff: 1.792}, 
UFFAtomType{name: "Li", atomic_symbol: "Li", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.336, theta: 1.91061, x: 3.021, d: 0.274, zeta: 13.969, z_eff: 2.703}, 
UFFAtomType{name: "Be3+2", atomic_symbol: "Be", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.074, theta: 1.91061, x: 4.499, d: 0.274, zeta: 13.969, z_eff: 2.703}, 
UFFAtomType{name: "B_3", atomic_symbol: "B", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.838, theta: 1.91061, x: 4.035, d: 0.274, zeta: 13.969, z_eff: 2.703}, 
UFFAtomType{name: "B_2", atomic_symbol: "B", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.828, theta: 1.60745, x: 4.035, d: 0.274, zeta: 13.969, z_eff: 2.703}, 
UFFAtomType{name: "C_3", atomic_symbol: "C", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.757, theta: 1.80118, x: 4.035, d: 0.274, zeta: 13.969, z_eff: 2.703}, 
UFFAtomType{name: "C_R", atomic_symbol: "C", bridging: false, aromatic: true, valency: 0, oxidation_state: 0, r: 0.729, theta: 1.91061, x: 4.035, d: 0.227, zeta: 14.866, z_eff: 2.348}, 
UFFAtomType{name: "C_2", atomic_symbol: "C", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.732, theta: 1.60919, x: 4.035, d: 0.044, zeta: 12.0, z_eff: 0.712}, 
UFFAtomType{name: "C_l", atomic_symbol: "C", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.706, theta: 2.09440, x: 3.947, d: 0.044, zeta: 12.0, z_eff: 0.712}, 
UFFAtomType{name: "N_3", atomic_symbol: "N", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.700, theta: 3.14159, x: 2.362, d: 0.056, zeta: 15.24, z_eff: 0.098}, 
UFFAtomType{name: "N_R", atomic_symbol: "N", bridging: false, aromatic: true, valency: 0, oxidation_state: 0, r: 0.699, theta: 1.57080, x: 2.451, d: 0.105, zeta: 12.73, z_eff: 1.912}, 
UFFAtomType{name: "N_2", atomic_symbol: "N", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.685, theta: 3.14159, x: 2.745, d: 0.105, zeta: 12.73, z_eff: 1.912}, 
UFFAtomType{name: "N_l", atomic_symbol: "N", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.656, theta: 1.91061, x: 4.083, d: 0.105, zeta: 12.73, z_eff: 1.912}, 
UFFAtomType{name: "O_3", atomic_symbol: "O", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.658, theta: 1.91061, x: 4.083, d: 0.025, zeta: 12.0, z_eff: 1.026}, 
UFFAtomType{name: "O_3_z", atomic_symbol: "O", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.528, theta: 2.09440, x: 3.851, d: 0.085, zeta: 12.0, z_eff: 1.565}, 
UFFAtomType{name: "O_R", atomic_symbol: "O", bridging: false, aromatic: true, valency: 0, oxidation_state: 0, r: 0.680, theta: 1.91061, x: 3.851, d: 0.180, zeta: 12.052, z_eff: 1.755}, 
UFFAtomType{name: "O_2", atomic_symbol: "O", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.634, theta: 2.09440, x: 3.851, d: 0.180, zeta: 12.052, z_eff: 1.755}, 
UFFAtomType{name: "O_1", atomic_symbol: "O", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 0.639, theta: 2.09440, x: 3.660, d: 0.069, zeta: 13.407, z_eff: 2.544}, 
UFFAtomType{name: "F_", atomic_symbol: "F", bridging: false, aromatic: false, valency: 1, oxidation_state: 0, r: 0.668, theta: 1.86227, x: 3.660, d: 0.069, zeta: 13.407, z_eff: 2.544}, 
UFFAtomType{name: "Ne4+4", atomic_symbol: "Ne", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 0.920, theta: 2.09440, x: 3.660, d: 0.069, zeta: 13.407, z_eff: 2.544}, 
UFFAtomType{name: "Na", atomic_symbol: "Na", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.539, theta: 1.94081, x: 3.660, d: 0.069, zeta: 13.407, z_eff: 2.544}, 
UFFAtomType{name: "A13", atomic_symbol: "A1", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.421, theta: 3.14159, x: 3.500, d: 0.060, zeta: 14.085, z_eff: 2.300}, 
UFFAtomType{name: "Si3", atomic_symbol: "Si", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.244, theta: 1.82404, x: 3.500, d: 0.060, zeta: 14.085, z_eff: 2.300}, 
UFFAtomType{name: "S_3+2", atomic_symbol: "S", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.117, theta: 2.54818, x: 3.500, d: 0.060, zeta: 14.085, z_eff: 2.300}, 
UFFAtomType{name: "S_3+4", atomic_symbol: "S", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.056, theta: 1.91986, x: 3.500, d: 0.060, zeta: 14.085, z_eff: 2.300}, 
UFFAtomType{name: "S_3+6", atomic_symbol: "S", bridging: false, aromatic: false, valency: 0, oxidation_state: 6, r: 1.064, theta: 2.09440, x: 3.500, d: 0.060, zeta: 14.085, z_eff: 2.300}, 
UFFAtomType{name: "S_R", atomic_symbol: "S", bridging: false, aromatic: true, valency: 0, oxidation_state: 0, r: 1.049, theta: 3.14159, x: 3.364, d: 0.050, zeta: 14.762, z_eff: 1.735}, 
UFFAtomType{name: "S_2", atomic_symbol: "S", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.027, theta: 3.14159, x: 3.243, d: 0.042, zeta: 15.440, z_eff: 0.194}, 
UFFAtomType{name: "Cl", atomic_symbol: "Cl", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.077, theta: 1.57080, x: 2.983, d: 0.030, zeta: 12.0, z_eff: 1.081}, 
UFFAtomType{name: "Ar4+4", atomic_symbol: "Ar", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 0.854, theta: 3.14159, x: 4.147, d: 0.305, zeta: 13.072, z_eff: 2.863}, 
UFFAtomType{name: "K_", atomic_symbol: "K", bridging: false, aromatic: false, valency: 1, oxidation_state: 0, r: 1.044, theta: 1.63712, x: 4.147, d: 0.305, zeta: 13.072, z_eff: 2.863}, 
UFFAtomType{name: "P_3+q", atomic_symbol: "P", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.032, theta: 1.91061, x: 4.147, d: 0.305, zeta: 13.072, z_eff: 2.863}, 
UFFAtomType{name: "Mg3+2", atomic_symbol: "Mg", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.953, theta: 1.91061, x: 4.295, d: 0.402, zeta: 12.175, z_eff: 2.323}, 
UFFAtomType{name: "P_3+3", atomic_symbol: "P", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.101, theta: 1.57080, x: 3.868, d: 0.185, zeta: 15.763, z_eff: 0.300}, 
UFFAtomType{name: "P_3+5", atomic_symbol: "P", bridging: false, aromatic: false, valency: 0, oxidation_state: 5, r: 1.056, theta: 3.14159, x: 3.812, d: 0.035, zeta: 12.0, z_eff: 1.165}, 
UFFAtomType{name: "Ca6+2", atomic_symbol: "Ca", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.761, theta: 1.57080, x: 3.399, d: 0.238, zeta: 12.0, z_eff: 2.141}, 
UFFAtomType{name: "Sc3+3", atomic_symbol: "Sc", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.513, theta: 1.91061, x: 3.144, d: 0.019, zeta: 12.0, z_eff: 2.592}, 
UFFAtomType{name: "TÍ3+4", atomic_symbol: "TÍ", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.412, theta: 1.91061, x: 3.295, d: 0.017, zeta: 12.0, z_eff: 2.659}, 
UFFAtomType{name: "Ti6+4", atomic_symbol: "Ti", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.412, theta: 1.57080, x: 3.175, d: 0.017, zeta: 12.0, z_eff: 2.659}, 
UFFAtomType{name: "V_3+5", atomic_symbol: "V", bridging: false, aromatic: false, valency: 0, oxidation_state: 5, r: 1.402, theta: 1.60745, x: 3.175, d: 0.016, zeta: 12.0, z_eff: 2.679}, 
UFFAtomType{name: "Cr6+3", atomic_symbol: "Cr", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.345, theta: 1.57080, x: 3.023, d: 0.015, zeta: 12.0, z_eff: 2.43}, 
UFFAtomType{name: "Mn6+2", atomic_symbol: "Mn", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.382, theta: 1.91061, x: 2.961, d: 0.005, zeta: 12.0, z_eff: 1.756}, 
UFFAtomType{name: "Fe3+2", atomic_symbol: "Fe", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.270, theta: 1.57080, x: 2.912, d: 0.124, zeta: 12.0, z_eff: 1.308}, 
UFFAtomType{name: "Fe6+2", atomic_symbol: "Fe", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.335, theta: 1.91061, x: 2.912, d: 0.415, zeta: 11.0, z_eff: 1.821}, 
UFFAtomType{name: "Co6+3", atomic_symbol: "Co", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.241, theta: 1.57080, x: 2.872, d: 0.379, zeta: 12.0, z_eff: 2.789}, 
UFFAtomType{name: "NÍ4+2", atomic_symbol: "NÍ", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.164, theta: 1.57080, x: 2.834, d: 0.309, zeta: 13.0, z_eff: 2.864}, 
UFFAtomType{name: "Cu3+1", atomic_symbol: "Cu", bridging: false, aromatic: false, valency: 0, oxidation_state: 1, r: 1.302, theta: 1.57080, x: 3.495, d: 0.291, zeta: 14.0, z_eff: 2.764}, 
UFFAtomType{name: "Zn3+2", atomic_symbol: "Zn", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.193, theta: 1.91061, x: 2.763, d: 0.251, zeta: 15.0, z_eff: 2.519}, 
UFFAtomType{name: "Ga3+3", atomic_symbol: "Ga", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.260, theta: 1.57080, x: 4.383, d: 0.015, zeta: 12.0, z_eff: 2.463}, 
UFFAtomType{name: "Ge3", atomic_symbol: "Ge", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.197, theta: 1.57080, x: 4.280, d: 0.013, zeta: 12.0, z_eff: 2.43}, 
UFFAtomType{name: "As3+3", atomic_symbol: "As", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.211, theta: 1.91061, x: 4.230, d: 0.013, zeta: 12.0, z_eff: 2.43}, 
UFFAtomType{name: "Se3+2", atomic_symbol: "Se", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.190, theta: 1.91061, x: 4.205, d: 0.013, zeta: 12.0, z_eff: 2.43}, 
UFFAtomType{name: "Br", atomic_symbol: "Br", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.192, theta: 1.91061, x: 4.189, d: 0.014, zeta: 12.0, z_eff: 2.43}, 
UFFAtomType{name: "Kr4+4", atomic_symbol: "Kr", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.147, theta: 1.91061, x: 4.141, d: 0.220, zeta: 16.0, z_eff: 0.452}, 
UFFAtomType{name: "Rb", atomic_symbol: "Rb", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 2.260, theta: 1.58127, x: 4.114, d: 0.04, zeta: 12.0, z_eff: 1.592}, 
UFFAtomType{name: "Sr6+2", atomic_symbol: "Sr", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 2.052, theta: 3.14159, x: 3.641, d: 0.235, zeta: 12.0, z_eff: 2.449}, 
UFFAtomType{name: "Y_3+3", atomic_symbol: "Y", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.698, theta: 1.57080, x: 3.345, d: 0.072, zeta: 12.0, z_eff: 3.257}, 
UFFAtomType{name: "Zr3+4", atomic_symbol: "Zr", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.564, theta: 3.14159, x: 3.124, d: 0.069, zeta: 12.0, z_eff: 3.667}, 
UFFAtomType{name: "Nb3+5", atomic_symbol: "Nb", bridging: false, aromatic: false, valency: 0, oxidation_state: 5, r: 1.473, theta: 1.57080, x: 3.165, d: 0.059, zeta: 12.0, z_eff: 3.618}, 
UFFAtomType{name: "M06+6", atomic_symbol: "M0", bridging: false, aromatic: false, valency: 0, oxidation_state: 6, r: 1.467, theta: 1.91061, x: 3.052, d: 0.056, zeta: 12.0, z_eff: 3.40}, 
UFFAtomType{name: "Mo3+6", atomic_symbol: "Mo", bridging: false, aromatic: false, valency: 0, oxidation_state: 6, r: 1.484, theta: 1.91061, x: 3.052, d: 0.056, zeta: 12.0, z_eff: 3.40}, 
UFFAtomType{name: "Tc6+5", atomic_symbol: "Tc", bridging: false, aromatic: false, valency: 0, oxidation_state: 5, r: 1.322, theta: 1.91061, x: 2.998, d: 0.048, zeta: 12.0, z_eff: 3.40}, 
UFFAtomType{name: "Xe4+4", atomic_symbol: "Xe", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.267, theta: 1.57080, x: 2.963, d: 0.007, zeta: 12.0, z_eff: 3.416}, 
UFFAtomType{name: "Ho6+3", atomic_symbol: "Ho", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.696, theta: 1.57080, x: 2.929, d: 0.007, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Er6+3", atomic_symbol: "Er", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.673, theta: 1.57080, x: 2.899, d: 0.006, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Ta3+5", atomic_symbol: "Ta", bridging: false, aromatic: false, valency: 0, oxidation_state: 5, r: 1.511, theta: 1.59872, x: 4.420, d: 0.228, zeta: 12.0, z_eff: 2.618}, 
UFFAtomType{name: "W_6+6", atomic_symbol: "W", bridging: false, aromatic: false, valency: 0, oxidation_state: 6, r: 1.392, theta: 1.57080, x: 4.404, d: 0.041, zeta: 12.0, z_eff: 3.271}, 
UFFAtomType{name: "W_3+4", atomic_symbol: "W", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.526, theta: 1.57080, x: 3.409, d: 0.056, zeta: 12.0, z_eff: 3.40}, 
UFFAtomType{name: "W_3+6", atomic_symbol: "W", bridging: false, aromatic: false, valency: 0, oxidation_state: 6, r: 1.380, theta: 1.57080, x: 3.391, d: 0.053, zeta: 12.0, z_eff: 3.508}, 
UFFAtomType{name: "Re6+5", atomic_symbol: "Re", bridging: false, aromatic: false, valency: 0, oxidation_state: 5, r: 1.372, theta: 1.91061, x: 3.170, d: 0.048, zeta: 12.0, z_eff: 3.21}, 
UFFAtomType{name: "Re3+7", atomic_symbol: "Re", bridging: false, aromatic: false, valency: 0, oxidation_state: 7, r: 1.314, theta: 1.57080, x: 3.069, d: 0.036, zeta: 12.0, z_eff: 1.956}, 
UFFAtomType{name: "Ru6+2", atomic_symbol: "Ru", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.478, theta: 1.91061, x: 3.069, d: 0.228, zeta: 12.0, z_eff: 1.65}, 
UFFAtomType{name: "Rh6+3", atomic_symbol: "Rh", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.332, theta: 1.91061, x: 3.069, d: 0.599, zeta: 11.0, z_eff: 2.07}, 
UFFAtomType{name: "Pd4+2", atomic_symbol: "Pd", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.338, theta: 1.57080, x: 2.954, d: 0.339, zeta: 12.0, z_eff: 2.961}, 
UFFAtomType{name: "Sb3+3", atomic_symbol: "Sb", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.407, theta: 1.91061, x: 2.954, d: 0.567, zeta: 13.0, z_eff: 2.704}, 
UFFAtomType{name: "Te3+2", atomic_symbol: "Te", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.386, theta: 1.57516, x: 4.470, d: 0.449, zeta: 14.0, z_eff: 2.882}, 
UFFAtomType{name: "I_", atomic_symbol: "I", bridging: false, aromatic: false, valency: 1, oxidation_state: 0, r: 1.382, theta: 3.14159, x: 4.50, d: 0.398, zeta: 15.0, z_eff: 2.65}, 
UFFAtomType{name: "Ag1+1", atomic_symbol: "Ag", bridging: false, aromatic: false, valency: 0, oxidation_state: 1, r: 1.386, theta: 3.14159, x: 3.148, d: 0.332, zeta: 12.0, z_eff: 0.556}, 
UFFAtomType{name: "Cd3+2", atomic_symbol: "Cd", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.403, theta: 1.91061, x: 2.848, d: 0.045, zeta: 12.0, z_eff: 1.573}, 
UFFAtomType{name: "In3+3", atomic_symbol: "In", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.459, theta: 1.91061, x: 4.463, d: 0.364, zeta: 12.0, z_eff: 2.727}, 
UFFAtomType{name: "Sn3", atomic_symbol: "Sn", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.398, theta: 1.91061, x: 4.392, d: 0.017, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Cs", atomic_symbol: "Cs", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 2.570, theta: 3.14159, x: 4.517, d: 0.013, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Ba6+2", atomic_symbol: "Ba", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 2.277, theta: 1.57080, x: 3.703, d: 0.010, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "La3+3", atomic_symbol: "La", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.943, theta: 1.91061, x: 3.522, d: 0.010, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Ce6+3", atomic_symbol: "Ce", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.841, theta: 1.57080, x: 3.556, d: 0.009, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Pr6+3", atomic_symbol: "Pr", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.823, theta: 1.57080, x: 3.606, d: 0.008, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Nd6+3", atomic_symbol: "Nd", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.816, theta: 1.57080, x: 3.575, d: 0.008, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Pm6+3", atomic_symbol: "Pm", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.801, theta: 1.57080, x: 3.547, d: 0.009, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Sm6+3", atomic_symbol: "Sm", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.780, theta: 1.57080, x: 3.520, d: 0.007, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Eu6+3", atomic_symbol: "Eu", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.771, theta: 1.57080, x: 3.493, d: 0.007, zeta: 12.0, z_eff: 3.30}, 
UFFAtomType{name: "Gd6+3", atomic_symbol: "Gd", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.735, theta: 1.57080, x: 3.368, d: 0.072, zeta: 12.0, z_eff: 3.921}, 
UFFAtomType{name: "Tb6+3", atomic_symbol: "Tb", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.732, theta: 1.57080, x: 3.451, d: 0.081, zeta: 12.0, z_eff: 4.075}, 
UFFAtomType{name: "Dy6+3", atomic_symbol: "Dy", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.710, theta: 1.57080, x: 3.428, d: 0.067, zeta: 12.0, z_eff: 3.70}, 
UFFAtomType{name: "Yb6+3", atomic_symbol: "Yb", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.637, theta: 1.57080, x: 3.355, d: 0.067, zeta: 12.0, z_eff: 3.70}, 
UFFAtomType{name: "Lu6+3", atomic_symbol: "Lu", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.671, theta: 1.57080, x: 3.640, d: 0.067, zeta: 12.0, z_eff: 3.70}, 
UFFAtomType{name: "Hf3+4", atomic_symbol: "Hf", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.611, theta: 1.91061, x: 3.141, d: 0.066, zeta: 12.0, z_eff: 3.70}, 
UFFAtomType{name: "Tm6+3", atomic_symbol: "Tm", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.660, theta: 1.57080, x: 3.374, d: 0.066, zeta: 12.0, z_eff: 3.70}, 
UFFAtomType{name: "Os6+6", atomic_symbol: "Os", bridging: false, aromatic: false, valency: 0, oxidation_state: 6, r: 1.372, theta: 1.57080, x: 3.120, d: 0.037, zeta: 12.0, z_eff: 3.70}, 
UFFAtomType{name: "Ir6+3", atomic_symbol: "Ir", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.371, theta: 1.57080, x: 2.840, d: 0.073, zeta: 12.0, z_eff: 3.731}, 
UFFAtomType{name: "Pt4+2", atomic_symbol: "Pt", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.364, theta: 1.57080, x: 2.754, d: 0.080, zeta: 12.0, z_eff: 3.382}, 
UFFAtomType{name: "Au4+3", atomic_symbol: "Au", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.262, theta: 1.57080, x: 3.293, d: 0.039, zeta: 12.0, z_eff: 2.625}, 
UFFAtomType{name: "Hgl+2", atomic_symbol: "Hg", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.340, theta: 3.14159, x: 2.705, d: 0.385, zeta: 12.0, z_eff: 1.75}, 
UFFAtomType{name: "T13+3", atomic_symbol: "T1", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.518, theta: 2.09440, x: 4.347, d: 0.680, zeta: 11.0, z_eff: 2.068}, 
UFFAtomType{name: "Pb3", atomic_symbol: "Pb", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.459, theta: 1.91061, x: 4.297, d: 0.663, zeta: 12.0, z_eff: 2.846}, 
UFFAtomType{name: "Ra6+2", atomic_symbol: "Ra", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 2.512, theta: 1.57080, x: 3.677, d: 0.518, zeta: 13.0, z_eff: 2.470}, 
UFFAtomType{name: "Ac6+3", atomic_symbol: "Ac", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.983, theta: 1.57080, x: 3.478, d: 0.325, zeta: 14.0, z_eff: 2.33}, 
UFFAtomType{name: "Th6+4", atomic_symbol: "Th", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.721, theta: 1.57080, x: 3.396, d: 0.284, zeta: 15.0, z_eff: 2.24}, 
UFFAtomType{name: "Pa6+4", atomic_symbol: "Pa", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.711, theta: 1.57080, x: 3.424, d: 0.248, zeta: 16.0, z_eff: 0.583}, 
UFFAtomType{name: "U_6+4", atomic_symbol: "U", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.684, theta: 1.57080, x: 3.395, d: 0.014, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Np6+4", atomic_symbol: "Np", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.666, theta: 1.57080, x: 3.424, d: 0.050, zeta: 12.0, z_eff: 1.847}, 
UFFAtomType{name: "Pu6+4", atomic_symbol: "Pu", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.657, theta: 1.57080, x: 3.424, d: 0.404, zeta: 12.0, z_eff: 2.92}, 
UFFAtomType{name: "Am6+4", atomic_symbol: "Am", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.660, theta: 1.57080, x: 3.381, d: 0.033, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Bi3+3", atomic_symbol: "Bi", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.512, theta: 1.57080, x: 4.370, d: 0.026, zeta: 12.0, z_eff: 4.202}, 
UFFAtomType{name: "Po3+2", atomic_symbol: "Po", bridging: false, aromatic: false, valency: 0, oxidation_state: 2, r: 1.50, theta: 1.57080, x: 4.709, d: 0.022, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "At", atomic_symbol: "At", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 1.545, theta: 3.14159, x: 4.750, d: 0.022, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Rn4+4", atomic_symbol: "Rn", bridging: false, aromatic: false, valency: 0, oxidation_state: 4, r: 1.420, theta: 1.57080, x: 4.765, d: 0.019, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Fr", atomic_symbol: "Fr", bridging: false, aromatic: false, valency: 0, oxidation_state: 0, r: 2.880, theta: 3.14159, x: 4.90, d: 0.016, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Cm6+3", atomic_symbol: "Cm", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.801, theta: 1.57080, x: 3.326, d: 0.013, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Bk6+3", atomic_symbol: "Bk", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.761, theta: 1.57080, x: 3.339, d: 0.013, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Cf6+3", atomic_symbol: "Cf", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.750, theta: 1.57080, x: 3.313, d: 0.013, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Es6+3", atomic_symbol: "Es", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.724, theta: 1.57080, x: 3.299, d: 0.012, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Fm6+3", atomic_symbol: "Fm", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.712, theta: 1.57080, x: 3.286, d: 0.012, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Md6+3", atomic_symbol: "Md", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.689, theta: 1.57080, x: 3.274, d: 0.011, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "No6+3", atomic_symbol: "No", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.679, theta: 1.57080, x: 3.248, d: 0.011, zeta: 12.0, z_eff: 3.90}, 
UFFAtomType{name: "Lw6+3", atomic_symbol: "Lw", bridging: false, aromatic: false, valency: 0, oxidation_state: 3, r: 1.698, theta: 1.57080, x: 3.236, d: 0.011, zeta: 12.0, z_eff: 3.90}, 
];
