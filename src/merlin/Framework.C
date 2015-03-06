#include <fstream>
#include <iostream>
#include <cstring>

#include <unistd.h>

using namespace std;
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "LatticeStructure.H"

#include "Vertex.H"
#include "Graph.H"

#include "FactorGraph.H"
#include "FactorManager.H"
#include "PotentialManager.H"
#include "MetaMove.H"
#include "MotifManager.H"
#include "BFGSWrapperData.H"
#include "BFGSWrapper.H"
#include "MetaLearner.H"

#include "Framework.H"

Framework::Framework()
{
	epsThreshold=-1;
}

Framework::~Framework()
{
}

//We will use getopt here
//The options are 
//-m modelname
//-o outputdir
//-e epsilon to control the number of standard deviations above random
//-k maxfactorsize
//-s number of samples for approximate information estimation
//-x k for which we approximate information
//-n cnt of the top candidate MBs to save

Error::ErrorCode
Framework::init(int argc, char** argv)
{
	int optret='-';
	opterr=1;
	int oldoptind=optind;
	int condCnt=1;
	while(optret=getopt(argc,argv,"m:o:k:d:u:v:l:p:q:r:s:t:c:h:")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return Error::UNKNOWN;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind ==2)
		{
			my_optarg=argv[oldoptind+1];	
		}
		else
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 'm':
			{
				char fName[256];
				sprintf(fName,"%s.model",my_optarg);
				Error::ErrorCode eCode=varManager.readVariables(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				evManager.setVariableManager(&varManager);
				sprintf(fName,"%s.data",my_optarg);
				eCode=evManager.loadEvidenceFromFile_Continuous(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				metaLearner.setGlobalEvidenceManager(&evManager);
				metaLearner.setVariableManager(&varManager);
				break;
			}
			case 'o':
			{
				metaLearner.setOutputDirName(my_optarg);
				break;
			}
			case 'k':
			{
				int aSize=atoi(my_optarg);
				metaLearner.setMaxFactorSize(1);
				metaLearner.setMaxFactorSize_Approx(aSize);
				break;
			}
			case 'd':
			{
				double convThreshold=atof(my_optarg);
				metaLearner.setConvergenceThreshold(convThreshold);
				break;
			}
			case 'u':
			{
				double lambda=atof(my_optarg);
				cout <<"Set Lambda" << lambda<< endl;
				metaLearner.setLambda(lambda);
				break;
			}
			case 'v':
			{
				cvCnt=atoi(my_optarg);
				break;
			}
			case 'l':
			{
				metaLearner.setRestrictedList(my_optarg);
				break;
			}
			case 'p':
			{
				metaLearner.setBeta1(atof(my_optarg));
				break;
			}
			case 'q':
			{
				metaLearner.setBeta_ChIP(atof(my_optarg));
				break;
			}
			case 'r':
			{
				metaLearner.setBeta_Motif(atof(my_optarg));
				break;
			}
			case 's':
			{
				metaLearner.setChIPGraph(my_optarg);
				break;
			}
			case 't':
			{
				metaLearner.setMotifGraph(my_optarg);
				break;
			}
			case 'c':
			{
				metaLearner.readModuleMembership(my_optarg);
				break;
			}
			case 'h':	
			{
				metaLearner.setClusteringThreshold(atof(my_optarg));
				break;
			}
			default:
			{
				cout <<"Unhandled option " << c  << endl;
				return Error::UNKNOWN;
			}
		}
		oldoptind=optind;
	}
	metaLearner.initPartitions(condCnt);
	return Error::SUCCESS;
}

int 
Framework::start()
{
	//metaLearner.start();
	metaLearner.doCrossValidation(cvCnt);
        return 0;
}


int
main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout <<"factorGraphInf " <<  endl
			<<"-m model" << endl
			<<"-o outputdir" << endl
			<< "-k maxfactorsize " << endl
			 << "-u lambda" << endl
			 << "-t convergence_threshold" << endl
			 << "-v cross_validation_cnt" << endl
			 << "-l restrictedfname" << endl
			 << "-p beta1" << endl
			 << "-q beta_chip" << endl
			 << "-r beta_motif" << endl
			 << "-s chipgraph" << endl
			 << "-t motifgraph" << endl
			<< "-c genelist"<< endl;

		return 0;
	}
	Framework fw;
	if(fw.init(argc,argv)!=Error::SUCCESS)
	{
		return 0;
	}
	fw.start();
	return 0;

}

