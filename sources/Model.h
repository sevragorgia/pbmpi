
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/



#include "RASCATGTRSBDPGammaPhyloProcess.h"

/*sevra
 * to reinclude the following add # before include
include "RASCATGTRFiniteGammaPhyloProcess.h"
include "RASCATFiniteGammaPhyloProcess.h"
include "RASCATSBDPGammaPhyloProcess.h"
include "AACodonMutSelSBDPPhyloProcess.h"
include "AACodonMutSelFinitePhyloProcess.h"
include "CodonMutSelFinitePhyloProcess.h"
include "CodonMutSelSBDPPhyloProcess.h"

#include "CodonSequenceAlignment.h"


*/

#include "Parallel.h"
#include <iostream>
#include <fstream>
#include <limits>

using namespace std;

MPI_Datatype Propagate_arg;



class Model	{

	public:

	PhyloProcess* process;
	string type;
	string name;
	int every;
	int until;
	int saveall;
	int incinit;

	Model(string datafile, string treefile, int modeltype, int nratecat, int mixturetype, int ncat, GeneticCodeType codetype, int suffstat, int fixncomp, int empmix, string mixtype, string rrtype, int iscodon, int fixtopo, int NSPR, int NNNI, int fixcodonprofile, int fixomega, int fixbl, int omegaprior, int kappaprior, int dirweightprior, double mintotweight, int dc, int inevery, int inuntil, int insaveall, int inincinit, int topoburnin, string inname, int myid, int nprocs)	{

		every = inevery;
		until = inuntil;
		name = inname;
		saveall = insaveall;
		incinit = inincinit;

		// 1 : CAT
		// 2 : CATGTR
		// 3 : MutSel

		// mixturetype
		// 0 : one for all
		// 1 : finite
		// 2 : dp
		// 3 : tdp : removed
		// 4 : sbdp
		// 5 : site specific
    
		/* sevra
     * original code with tons of options...
		// CAT
		if (modeltype == 1)	{
			if (myid == 0) {
				// cerr << "cat model\n";
			}
			if (mixturetype == 1)	{
				type = "CATFINITE";
				process = new RASCATFiniteGammaPhyloProcess(datafile,treefile,nratecat,ncat,fixncomp,empmix,mixtype,dirweightprior,fixtopo,NSPR,NNNI,dc,myid,nprocs); 
			}
			else	{
				type = "CATSBDP";
				process = new RASCATSBDPGammaPhyloProcess(datafile,treefile,nratecat,iscodon,codetype,fixtopo,NSPR,NNNI,kappaprior,dirweightprior,mintotweight,dc,incinit,myid,nprocs); 
			}
		}

		// CATGTR
		else if (modeltype == 2)	{
			if (myid == 0) {
				// cerr << "catgtr model\n";
			}
			if (mixturetype == 1)	{
				if (suffstat)	{
					type = "CATGTRFINITE";
					process = new RASCATGTRFiniteGammaPhyloProcess(datafile,treefile,nratecat,ncat,fixncomp,empmix,mixtype,rrtype,dirweightprior,fixtopo,NSPR,NNNI,dc,myid,nprocs); 
				}
				else	{
					cerr << "gpss deprecated\n";
					exit(1);
				}
			}
			else if (mixturetype == 2)	{
				cerr << "simple dp deprecated\n";
				exit(1);
			}
			else if (mixturetype == 3)	{
				if (suffstat)	{
					type = "CATGTRSBDP";
					process = new RASCATGTRSBDPGammaPhyloProcess(datafile,treefile,nratecat,iscodon,codetype,rrtype,fixtopo,NSPR,NNNI,kappaprior,dirweightprior,mintotweight,dc,incinit,myid,nprocs); 
				}
				else	{
					cerr << "gpss deprecated\n";
					exit(1);
				}
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized\n";
				exit(1);
			}
		}

		// AAMUTSEL
		else if (modeltype == 3)	{
			cerr << "deprecated.\n";
			exit(1);
		}
		

		// CodonMutSel
		else if (modeltype == 4)	{
			if (mixturetype == 1)	{
				type = "CODONMUTSELFINITE";
				process = new CodonMutSelFinitePhyloProcess(datafile,treefile,codetype,ncat,fixncomp,empmix,mixtype,fixtopo,fixbl,NSPR,NNNI,dirweightprior,dc,myid,nprocs);
			}
			else if (mixturetype == 3)	{
				type = "CODONMUTSELSBDP";
				process = new CodonMutSelSBDPPhyloProcess(datafile,treefile,codetype,fixtopo,fixbl,NSPR,NNNI,kappaprior,mintotweight,dc,myid,nprocs);
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized or not yet implemented.\n";
				exit(1);
			}
		}

		// AACodonMutSel
		else	{
			if (mixturetype == 1)	{
				type = "AACODONMUTSELFINITE";
				process = new AACodonMutSelFinitePhyloProcess(datafile,treefile,codetype,ncat,fixncomp,empmix,mixtype,fixtopo,fixbl,NSPR,NNNI,fixcodonprofile,fixomega,omegaprior,dirweightprior,dc,myid,nprocs);
			}
			else if (mixturetype == 3)	{
				type = "AACODONMUTSELSBDP";
				process = new AACodonMutSelSBDPPhyloProcess(datafile,treefile,codetype,fixtopo,fixbl,NSPR,NNNI,fixcodonprofile,fixomega,omegaprior,kappaprior,dirweightprior,mintotweight,dc,myid,nprocs);
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized or not yet implemented.\n";
				exit(1);
			}
			
		}
		*/
		
    // CATGTR only code
		if (modeltype == 2)	{
			if (myid == 0) {
				// cerr << "catgtr model\n";
			}
			
			if (mixturetype == 3)	{
				if (suffstat)	{
					type = "CATGTRSBDP";
					process = new RASCATGTRSBDPGammaPhyloProcess(datafile,treefile,nratecat,iscodon,codetype,rrtype,fixtopo,NSPR,NNNI,kappaprior,dirweightprior,mintotweight,dc,incinit,myid,nprocs); 
				}
				else	{
					cerr << "gpss deprecated\n";
					exit(1);
				}
			}else{//sevra
        
        cerr << "only mixture possible is GTRSBDP \n";
        exit(1);
        
      }
    } else {
      
        cerr << "only CATGTR is available \n";
        exit(1);
      
    }
    
		process->SetFixBL(fixbl);
		if (fixbl && (! myid))	{
			cerr << "set lengths from tree file\n";
			process->SetLengthsFromNames();
		}
		process->SetTopoBurnin(topoburnin);
	}

	Model(string inname, int myid, int nprocs)	{

		name = inname;

		ifstream is((name + ".param").c_str());
		if (! is)	{
			cerr << "error: cannot open " << name << ".param\n";
			exit(1);
		}

		is >> type;
		int size;
		is >> every >> until >> size;
		is >> saveall;
		
    /* sevra
     * original code with lots of options
		if (type == "CATSBDP")	{
			process = new RASCATSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATFINITE")	{
			process = new RASCATFiniteGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATGTRSBDP")	{
			process = new RASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATGTRFINITE")	{
			process = new RASCATGTRFiniteGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "AACODONMUTSELFINITE")	{
			process = new AACodonMutSelFinitePhyloProcess(is,myid,nprocs);
		}
		else if (type == "CODONMUTSELFINITE")	{
			process = new CodonMutSelFinitePhyloProcess(is,myid,nprocs);
		}
		else if (type == "CODONMUTSELSBDP")	{
			process = new CodonMutSelSBDPPhyloProcess(is,myid,nprocs);
		}
		else if (type == "AACODONMUTSELSBDP")	{
			process = new AACodonMutSelSBDPPhyloProcess(is,myid,nprocs);
		}
		else	{
			cerr << "error, does not recognize model type : " << type << '\n';
			exit(1);
		}
    */
    
    if (type == "CATGTRSBDP")	{//we only want CATGTRSBDP
			process = new RASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
    
		process->SetSize(size);
	}

	void ToStream(ostream& os, bool header)	{
		stringstream ss;
		if (header)	{
			ss << type << '\n';
			ss << every << '\t' << until << '\t' << GetSize() << '\n';
			ss << saveall << '\n';
			process->ToStreamHeader(ss);
		}
		process->ToStream(ss);
		os << ss.str();
	}

	~Model()	{
		// things here !!
	}

	void WaitLoop()	{
		process->WaitLoop();
	}

	double Move(double tuning, int nrep)	{
		double total = 0;
    //sevra
    cout << "in Move 1 " << nrep << "\n";
		for (int rep=0; rep<nrep; rep++)	{
			total += process->Move(tuning);
		}
		cout << total/nrep << "\n";
		return total / nrep;
	}

	int RunningStatus()	{
		ifstream ris((name + ".run").c_str());
		int i;
		ris >> i;
		ris.close();
		return i;
	}

	int GetSize()	{
		return process->GetSize();
	}
	
	void IncSize()	{
		process->IncSize();
	}

	void Run(int burnin)	{
		if (burnin != 0)	{
			if (GetSize() < burnin)	{
				process->SetBurnin(true);
			}
		}
		ofstream ros((name + ".run").c_str()); stringstream buf;
		buf << 1 << '\n';
		ros << buf.str();
		ros.close();
	
		while (RunningStatus() && ((until == -1) || (GetSize() < until)))	{
			if (GetSize() >= burnin)	{
				process->SetBurnin(false);
			}

			//sevra
			cout << "in Run\n";
			Move(1,every);
			
			process->IncSize();

			ofstream os((name + ".treelist").c_str(), ios_base::app);
			TreeTrace(os);
			os.close();

			ofstream tos((name + ".trace").c_str(), ios_base::app);
			Trace(tos);
			tos.close();

			ofstream mos((name + ".monitor").c_str());
			Monitor(mos);
			mos.close();

			ofstream pos((name + ".param").c_str());
			pos.precision(numeric_limits<double>::digits10);
			ToStream(pos,true);
			pos.close();

			if (saveall)	{
				ofstream cos((name + ".chain").c_str(),ios_base::app);
				cos.precision(numeric_limits<double>::digits10);
				ToStream(cos,false);
				cos.close();
			}

		}	
		cerr << name << ": stopping after " << GetSize() << " points.\n";
		cerr << '\n';
	}

	NewickTree* GetTree() {return process->GetLengthTree();}

	void TraceHeader(ostream& os)	{
		stringstream ss;
		process->TraceHeader(ss);
		os << ss.str();
	}

	void Trace(ostream& os)	{
		stringstream ss;
		process->Trace(ss);
		os << ss.str();
	}
	
	void TreeTrace(ostream& os)	{
		stringstream ss;
		process->SetNamesFromLengths();
		process->RenormalizeBranchLengths();
		GetTree()->ToStream(ss);
		process->DenormalizeBranchLengths();
		os << ss.str();
	}
	
	void Monitor(ostream& os)	{
		stringstream ss;
		process->Monitor(ss);
		os << ss.str();
	}

	void ReadPB(int argc, char* argv[])	{
		process->ReadPB(argc,argv);
	}
};
