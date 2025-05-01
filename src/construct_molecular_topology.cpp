#include<iostream> 
#include<string> 
#include<armadillo>
#include<cmath>
#include "math_functions.h"

/*
STRUCTS AND CLASSES USED TO CONSTRUCT THE MOLECULAR TOPOLOGY:
*/
struct Bond {
	double distance; 
	double atom1;
	double atom2; 
	// add other constants later e.g. force constant etc. 
};


class Molecule {
public:
        // Attributes
        int nr_atoms;
        arma::vec atom_vec;
	arma::vec enumerated_atom_vec; 
        arma::mat xyz_mat;
	std::vector<std::pair<double, arma::vec>> neighbor_vec;
	std::vector<Bond> bond_vec; 
	arma::vec bond_angles; 
	arma::mat out_of_plane_angles; 
	arma::mat torsiona_angles; 

        // Methods
        void makeNeighborList () {
                double rH = 31 * 0.01;  					   			 // pm * 0.01 Å/pm, covalent radius H
                double rC = 76 * 0.01;  					  			 // pm * 0.01 Å/pm, covalent radius sp^3 C

                for (int i = 0; i < nr_atoms; i++) {
                        std::vector<std::pair<double, arma::vec>> neighborlst_i;
                        arma::vec neighbors; 

                        for (int j = 0; j < nr_atoms; j++) {  // two first C-H bonds, last C-C bonds
                                if (euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rH + rC) && enumerated_atom_vec(j) != enumerated_atom_vec(i) && atom_vec(i) == 6 && atom_vec(j) == 1                                     || 																									    euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rC + rC) && enumerated_atom_vec(j) != enumerated_atom_vec(i) && atom_vec(i) == 6 && atom_vec(j) == 6                                           ||  																							            euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rH + rC) && enumerated_atom_vec(j) != enumerated_atom_vec(i) && atom_vec(i) == 1 && atom_vec(j) == 6)                                 {                       			
				    neighbors.insert_rows(neighbors.n_rows, arma::rowvec({enumerated_atom_vec(j)}));  
                                 }
                        }
                        std::pair<double, arma::vec> entry = {enumerated_atom_vec(i), neighbors};
                        neighbor_vec.push_back(entry);
                }
        }


	bool bondExists(double a1, double a2, const std::vector<Bond> bond_vec) {
        for (const Bond b : bond_vec) {
                if ((b.atom1 == a1 && b.atom2 == a2) || (b.atom1 == a2 && b.atom2 == a1)) {
                        return true;
                }
        }
        return false;
	}


	void makeChemicalBonds () {   // Want to loop over all the carbons and create bond from each carbon to its neighbour (but not doulbe count C-C bonds!)
		for (int i = 0; i < nr_atoms; i++) {
			for (int j : neighbor_vec[i].second) {
				if (i != j && neighbor_vec[i].second.n_elem > 1) { 					// Base the connection on the C-atoms 
					int atom1 = i;
					int atom2 = j;

					if (!bondExists(atom1, atom2, bond_vec)) {
						Bond bond; 
						bond.distance = euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t());
						bond.atom1 = enumerated_atom_vec(i);
						bond.atom2 = enumerated_atom_vec(j);
						bond_vec.push_back(bond);
					}
				}
                       	}
        	}
	}



	void getBondAngles () {
		for (int i = 0; i < bond_vec.size(); i++) {
			for (int j = 0; j < i; j++) {
				int center_atom; 
				if (bond_vec[i].atom1 == bond_vec[j].atom1) center_atom = bond_vec[i].atom1; 
				else if (bond_vec[i].atom2 == bond_vec[j].atom2) center_atom = bond_vec[i].atom2;
				else if (bond_vec[i].atom1 == bond_vec[j].atom2) center_atom = bond_vec[i].atom1;
				else if (bond_vec[i].atom2 == bond_vec[j].atom1) center_atom = bond_vec[i].atom2;
				else continue;  
					
					double external_atom1;
					double external_atom2;
					double angle;
		
					if (center_atom == bond_vec[i].atom1) {
						external_atom1 = bond_vec[i].atom2;
					} else {
						external_atom1 = bond_vec[i].atom1;
					}

					if (center_atom == bond_vec[j].atom1) {
						external_atom2 = bond_vec[j].atom2;
					} else {
						external_atom2 = bond_vec[j].atom1;
					}

					angle = acos(arma::dot(getUnitVector(xyz_mat.row(center_atom).t(), xyz_mat.row(external_atom1).t()), getUnitVector(xyz_mat.row(center_atom).t(), xyz_mat.row(external_atom2).t()))) * (180.0/acos(-1.0)); 
	                		//bond_angles.insert_rows(bond_angles.n_rows, angle);
								
					std::cout << external_atom1 << " " << center_atom << " " << external_atom2 << std::endl;
					std::cout << angle << std::endl;
				}
			}
		}	
};


/*
FUNCTIONS USED TO CONSTRUCT THE MOLECULAR TOPOLOGY:
*/

// Function used to read in atomic positions from a .txt file to an armadillo matrix. 
arma::mat readFile(std::string filename) {
	std::ifstream file (filename);
		if (file.is_open()) {
			std::string first_line;
        		std::getline(file, first_line);  // read and discard the first line
			
			int nr_atoms = std::stoi(first_line);
			arma::mat input_data_matrix = arma::mat(nr_atoms, 4); 

			for (int i = 0; i < nr_atoms; i++) {
				for (int j = 0; j < 4; j++) {
					file >> input_data_matrix(i,j);
				}
			}
			file.close(); 
			return input_data_matrix; 
		} else {
			std::cerr << "Error: Failed to open file." << std::endl;
	}
}


// Function used to enumerate an arma::vec object. 
arma::vec enumerateVector(arma::vec vector) {
	for (int i = 0; i < vector.n_elem; i++) {
		vector(i) = i;
	}
	return vector; 
}



// Main function that executes the program. 
int main(int argc, char** argv) {
	// Read data: 
	std::string filepath;

	if (argc < 2) {
		std::cerr << "Error: Format for input filepath must be:  " << argv[0] << " <file_path>\n";
		return 1; 
	} else {
		filepath = argv[1];
	}

	arma::mat input_mat = readFile(filepath);

	// Construct molecule: 
	Molecule mol;
	mol.nr_atoms = input_mat.n_rows; 
	mol.atom_vec = input_mat.col(0);
	mol.enumerated_atom_vec = enumerateVector(input_mat.col(0));
	mol.xyz_mat  = input_mat.cols(1,3); 
	mol.makeNeighborList(); 
        mol.makeChemicalBonds();
	mol.getBondAngles();

	
	std::cout << "mol.nr_atoms" << std::endl;
	std::cout << mol.nr_atoms << std::endl;
	std::cout << "mol.atom_vec" << std::endl;
	std::cout << mol.atom_vec << std::endl;
	std::cout << "xyz_mat"      << std::endl; 
	std::cout << mol.xyz_mat << std::endl;
	std::cout << "mol.neighbor_vec" << std::endl;

for (size_t i = 0; i < mol.neighbor_vec.size(); ++i) {
    const auto& entry = mol.neighbor_vec[i];
    std::cout << "Atom: " << entry.first << " has neighbors: ";

    for (size_t j = 0; j < entry.second.n_elem; ++j) {
        std::cout << entry.second(j);
        if (j != entry.second.n_elem - 1)
            std::cout << ", ";
    }

    std::cout << std::endl;
}

for (const auto& bond : mol.bond_vec) {
    std::cout << "Bond between atom1: " << bond.atom1 
              << " and atom2: " << bond.atom2 
              << " with distance: " << bond.distance << std::endl;
}

	std::cout << "bond angles" << std::endl;
	mol.bond_angles.print(); 

return 0;
}
