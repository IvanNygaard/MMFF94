// Function used to compute the distance between two atoms at positions vec1 and vec2: 
double euclideanDistance(arma::vec vec1, arma::vec vec2) {
	double distance; 
	distance = std::pow(std::pow(vec1(0) - vec2(0), 2) + std::pow(vec1(1) - vec2(1), 2) + std::pow(vec1(2) - vec2(2), 2), 0.5);
	
	return distance; 
}


// Function used to get the unit vector in the vec1 -> vec2 direction: 
arma::vec getUnitVector(arma::vec vec1, arma::vec vec2) { 
	arma::vec unitvector;
	double x_comp;
	double y_comp;
	double z_comp; 
	double distance;

	distance = euclideanDistance(vec1, vec2); 

	x_comp = -(vec1(0) - vec2(0)) * std::pow(distance, -1); 
	y_comp = -(vec1(1) - vec2(1)) * std::pow(distance, -1); 
	z_comp = -(vec1(2) - vec2(2)) * std::pow(distance, -1); 

	unitvector = {x_comp, y_comp, z_comp}; 

	return unitvector; 
}


// Function used to compute the factorial of an integer n.
int factorial(int n) {
        // Compute the factorial for an integer n >= 0.
        int res = 1;
        int k   = 1;

        while (k <= n) {
                res = res * k;
                k   = k + 1;
        }
        return res;
}


// Function used to compute the binomial coefficient (nCk): 
int binomialCoefficient(int n, int k) {
        return (factorial(n))/(factorial(k) * factorial(n-k));
}
