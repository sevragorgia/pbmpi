
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "BranchProcess.h"
#include "Random.h"

#include <sstream>
#include <iomanip>

string BranchProcess::GetNodeName(const Link* link)  const	{
	if (link->isLeaf())	{
		return link->GetNode()->GetName();
	}
	else	{
		return "";
	}
}

string BranchProcess::GetBranchName(const Link* link) const	{
	if (link->isRoot())	{
		return "";
	}
	else	{
		ostringstream s;
		s << GetLength(link->GetBranch());
		return s.str();
	}
}	

void BranchProcess::Backup()	{
	for (int i=0; i<GetNbranch(); i++)	{
		bkarray[i] = blarray[i];
	}
}

void BranchProcess::Restore()	{
	for (int i=0; i<GetNbranch(); i++)	{
		blarray[i] = bkarray[i];
	}
}

void BranchProcess::Restore(const Branch* branch)	{
	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	int index = branch->GetIndex();
	blarray[index] = bkarray[index];
}

void BranchProcess::MoveBranch(const Branch* branch, double m)	{

	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	
	int index = branch->GetIndex();
	
  //sevra: this is the method that is called to update branch lengths
	cout << "proposing bl move for branch " << index << ", m=" << m << ", " << exp(m) << "\n"; 
	
  bkarray[index] = blarray[index];
	blarray[index] *= exp(m);
}

double BranchProcess::ProposeMove(const Branch* branch, double tuning)	{

	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	
	//sevra: I think branches are updated here.
	cout << "proposing bl move: ";
  
	int index = branch->GetIndex();
	bkarray[index] = blarray[index];
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
  
  //sevra
  cout << "m,  " << exp(m) << "\n";
  
	blarray[index] *= exp(m);
	return m;
}

int BranchProcess::MoveAllBranches(double e)	{
	return RecursiveMoveAllBranches(GetRoot(),e);
}

int BranchProcess::RecursiveMoveAllBranches(const Link* from, double e)	{

	int n = 0;
	if (!from->isRoot())	{
		int index = from->GetBranch()->GetIndex();
		blarray[index] *= e;
		n++;
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		n += RecursiveMoveAllBranches(link->Out(),e);
	}
	return n;
}

double BranchProcess::RecursiveLogLengthPrior(const Link* from)	{

	double total = 0;
	if (! from->isRoot())	{
		total += LogBranchLengthPrior(from->GetBranch());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += RecursiveLogLengthPrior(link->Out());
	}
	return total;
}

void BranchProcess::RecursiveSampleLength(const Link* from)	{

	if (! from->isRoot())	{
		SampleLength(from->GetBranch());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSampleLength(link->Out());
	}
}

void BranchProcess::RecursiveNormalizeBranchLengths(const Link* from, double factor)	{

	if (! from->isRoot())	{
		SetLength(from->GetBranch(),GetLength(from->GetBranch())*factor);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveNormalizeBranchLengths(link->Out(),factor);
	}
}

double BranchProcess::RecursiveTotalLength(const Link* from)	{
	double total = 0;
	if (! from->isRoot())	{
		total += GetLength(from->GetBranch());
	}
	else	{
		if (GetLength(from->GetBranch()))	{
			cerr << "error: non null branch length for root\n";
			cerr << GetLength(from->GetBranch()) << '\n';
			exit(1);
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += RecursiveTotalLength(link->Out());
	}
	return total;
}

/*sevra: this set the length of the branches?*/
void BranchProcess::RecursiveSetLengthsFromNames(const Link* from)	{
	if (! from->isRoot())	{
    
    /*sevra: this was the original code
    
		double l = atof(from->GetBranch()->GetName().c_str());
		
		I wonder whether should be rather be:
    double l = atof(from->GetNode()->GetName().c_str());
     
     also, apparently atof will only result in a number if the text passed as an argument is a number.
     As a result, all branches get the same length 1e-10 after this method is called.
    
    I guess this method set the branch lengths when importing a tree from a file.
    
     */
    
    double l = atof(from->GetBranch()->GetName().c_str());
    
    cout << "attempting to set branch length: " << l << ", for branch " << from->GetBranch()->GetName() << "\n";
    
    
    if (l > 0){
      cout << "+l\n";
    }
    
		if (l < 0)	{
      cout << "-l\n";
			cerr << "error in BranchProcess::SetLengthsFromFile : negative branch length: " << l << '\n';
			exit(1);
		}
		if (l == 0)	{
      cout << "0=l\n";
      //sevra. All branches are set to this value for some strange reason... maybe change it to be sampled from a distribution?
			l = 1e-10;
			// cerr << "warning: null branch length\n";
		}
		/*
		if (l <= 0)	{
			// l = 0.1;
			cerr << "error in BranchProcess::SetLengthsFromFile : " << l << '\n';
			exit(1);
		}
		*/
    
    //cout << "setting branch " << from->GetBranch()->GetIndex() << ", " << l << "\n";
    
		SetLength(from->GetBranch(),l);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetLengthsFromNames(link->Out());
	}
}


void BranchProcess::RecursiveSetNamesFromLengths(const Link* from)	{

	if (! from->isRoot())	{
		ostringstream s;
		s.precision(15);
		s << GetLength(from->GetBranch());
		from->GetBranch()->SetName(s.str());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetNamesFromLengths(link->Out());
	}
}

double BranchProcess::LengthSuffStatLogProb()	{

	double total = 0;
	for (int i=1; i<GetNbranch(); i++)	{
		total += GetBranchLengthSuffStatCount(i) * log(blarray[i]);
		total -= GetBranchLengthSuffStatBeta(i) * blarray[i];
	}
	return total;
}



