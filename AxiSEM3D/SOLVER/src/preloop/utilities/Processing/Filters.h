// create filters 

#pragma once 

#include "eigenc.h"
#include <cmath>

class Filters{
	
public:
	
	static void logGabor(RMatXX &filters, const RColX &freq, Real center, Real sigma, int indFilt) {
			

			filters(indFilt,0) = 0.;
			
			for (int i = 1; i < freq.size(); i++) {
				
				filters(indFilt,i) = exp( (-pow(log(freq(i)/center),2.) )/(2*pow(log(sigma),2.)) );
				
			}
			
		
		
	}
	
	
	static void butterLowpass(RMatXX &filters, const RColX &freq, Real cutoff, int order, int indFilt) {
			
			
			for (int i = 0; i < freq.size(); i++) {
				
				filters(indFilt,i) = 1 / sqrt(1 + pow(freq(i)/cutoff, 2. * order));
				
			}
			
		
		
	}
	
};