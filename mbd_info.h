#ifndef MBD_INFO_H
#define MBD_INFO_H

#include "TH2Poly.h"

double a = 28.0/17.0;
double r = a*cos(30*3.141592/180);
double width = 28.0/2;
double height = 10*r;

double mbd_coords[64][2] = {
    {-12.34447187692347, 4.13084107008881  }, 
    {-12.416453509420549, 1.364486440539988}, 
    {-10.041133773551367, 8.392523393616344}, 
    {-9.933161324805749, 5.663551795985811 }, 
    {-9.861184634744298, 2.6728989749272642}, 
    {-7.521850772688031, 9.850468055676771 }, 
    {-7.4858599564394925, 6.934580015005359}, 
    {-7.4858599564394925, 4.242991449293115}, 
    {-7.413883266378041, 1.2897203767034124}, 
    {-5.002572714260323, 11.345794466206044}, 
    {-5.002572714260323, 8.467289457452917 }, 
    {-4.930591081763241, 5.700936111353542 }, 
    {-2.483289713396985, 12.76635481289874 }, 
    {-2.4472988971484426, 9.850468055676771}, 
    {-2.4113130233355324, 7.121495174596795}, 
    {-0.03598834503072368, 11.4205605300426},    
    {-0.03598834503072368, 8.46728945745291},
    {2.4472988971484426, 12.803737844817029}, 
    {2.5192805296455276, 9.962617151431633 }, 
    {2.4113130233355378, 7.08411214267851  }, 
    {4.894600265514704, 11.308411434287756 }, 
    {4.930591081763243, 8.467289457452917  }, 
    {4.930591081763243, 5.700936111353542  }, 
    {7.485859956439487, 9.925234119513348  }, 
    {7.485859956439487, 6.971961763474205  }, 
    {7.377892450129497, 4.280373197761957  }, 
    {7.413883266378047, 1.4018694724582765 }, 
    {9.933161324805745, 8.467289457452917  }, 
    {9.933161324805745, 5.663551795985811  }, 
    {9.861184634744301, 2.7850467872326874 }, 
    {12.416453509420549, 4.168225385456537 }, 
    {12.38046269317201, 1.364486440539988  }, 
    {12.38046269317201, -4.2429888823942274}, 
    {12.344471876923475, -1.401868189008832}, 
    {9.861184634744301, -8.542054237840048 }, 
    {9.897170508557206, -5.663549229086925 }, 
    {9.969152141054291, -2.8224285357015297}, 
    {7.341901633880962, -9.999998899900477 }, 
    {7.305910817632423, -7.084109575779622 }, 
    {7.449869140190952, -4.3551392615985325}, 
    {7.413883266378047, -1.5514003166819812}, 
    {4.930591081763243, -11.383174931225444}, 
    {4.930591081763243, -8.542054237840048 }, 
    {4.930591081763243, -5.663549229086925 }, 
    {2.4832897133969816, -12.84111959328587},    
    {2.4832897133969816, -9.962614584532748}, 
    {2.4832897133969816, -7.084109575779622}, 
    {-0.03598834503072368, -11.383174931225},    
    {-0.03598834503072368, -8.5420542378400},   
    {-2.5192805296455223, -12.8411195932858},   
    {-2.5192805296455223, -9.96261458453274},    
    {-2.5192805296455223, -7.08410957577962},    
    {-4.96658189801178, -11.383174931225444}, 
    {-4.96658189801178, -8.504672489371206 }, 
    {-4.894600265514699, -5.700933544454655}, 
    {-7.449874082626583, -9.925232836063902}, 
    {-7.449874082626583, -7.196259954983927}, 
    {-7.449874082626583, -4.242988882394227},    
    {-7.449874082626583, -1.439252504376561},    
    {-9.969152141054288, -8.542054237840048}, 
    {-10.041133773551367, -5.70093354445465},    
    {-10.041133773551367, -2.78504678723268},    
    {-12.488430199482, -4.2429888823942274 }, 
    {-12.488430199482, -1.401868189008832  }
};

int mbd_ring_index[64] = 
    {2,2,2,1,1,2,1,0,
     0,2,1,0,2,1,0,1,
     0,2,1,0,2,1,0,2,
     1,0,0,2,1,1,2,2,
     2,2,2,1,1,2,1,0,
     0,2,1,0,2,1,0,1,
     0,2,1,0,2,1,0,2,
     1,0,0,2,1,1,2,2};

float gaincorr[128] = { 0 };
float tq_t0_offsets[128] = { 0 };

void get_mbd_corrections() {
    std::ifstream gainfile( "gainfile.calib" );

    int ch;
    float integ, integerr;
    float peak, peakerr;
    float mbd_width, widtherr;
    float chi2ndf;
    while ( gainfile >> ch >> integ >> peak >> mbd_width >> integerr >> peakerr >> widtherr >> chi2ndf )
    {
        gaincorr[ch] = 1.0/peak;
    }

    gainfile.close();

    std::ifstream tfile( "/sphenix/user/samfred/commissioning/macros/calib/bbc_tq_t0.calib" );

    int pmtnum;
    float meanerr;
    float sigma;
    float sigmaerr;
    for (int ipmt = 0; ipmt < 128; ipmt++) {
        tfile >> pmtnum >> tq_t0_offsets[ipmt] >> meanerr >> sigma >> sigmaerr;
    }
}

void mbd_honeycomb(TH2Poly * hist, Double_t xstart, Double_t ystart, Double_t a,Int_t k, Int_t s) {

    // a is the side length of the hexagon
    
    Double_t numberOfHexagonsInTheRow;
    Double_t x[6], y[6];
    Double_t xloop, yloop, ytemp;
    xloop = xstart; yloop = ystart + a*TMath::Sqrt(3)/2.0;
    for (int sCounter = 0; sCounter < s; sCounter++) {

        ytemp = yloop; // Resets the temp variable

        // Determine the number of hexagons in that row
        if(sCounter%2 == 0){numberOfHexagonsInTheRow = k;}
        else{numberOfHexagonsInTheRow = k - 1;}

        for (int kCounter = 0; kCounter <  numberOfHexagonsInTheRow; kCounter++) {

            // Go around the hexagon
            x[0] = xloop;
            y[0] = ytemp;
            x[1] = x[0] + a/2.0;
            y[1] = y[0] + a*TMath::Sqrt(3)/2.0;
            x[2] = x[1] + a;
            y[2] = y[1];
            x[3] = x[2] + a/2.0;
            y[3] = y[0];
            x[4] = x[2];
            y[4] = y[3] - a*TMath::Sqrt(3)/2.0;
            x[5] = x[1];
            y[5] = y[4];

            hist->AddBin(6, x, y);

            // Go up
            ytemp += a*TMath::Sqrt(3);
        }

        // Increment the starting position
        if (sCounter%2 == 0) yloop += a*TMath::Sqrt(3)/2.0;
        else                 yloop -= a*TMath::Sqrt(3)/2.0;
        xloop += 1.5*a;
    }
}

map<int,vector<bool>> readmap(ifstream& inmap) {
    map<int,vector<bool>> map = {};
    if (inmap.is_open()) {
        string line;
        while (getline(inmap,line)) {
            istringstream iss(line);
            int count = 0;
            int runnum = 0;
            vector<bool> evts = {};
            do {
                string word;
                iss >> word;
                if (count > 0 && isdigit(*word.c_str())) evts.push_back(stoi(word));
                else if (isdigit(*word.c_str())) runnum = stoi(word);
                count++;
            } while(iss);
            map[runnum] = evts;
        }
    }
    return map;
}

map<int,vector<bool>> good_evts5 = {};
map<int,vector<bool>> good_evts30 = {};

void mbd_init() {
    get_mbd_corrections();

    //ifstream inmap5("event_map_5cm.txt");
    //ifstream inmap30("event_map_30cm.txt");
    //good_evts5 = readmap(inmap5);
    //good_evts30 = readmap(inmap30);

    return;
} 

#endif
