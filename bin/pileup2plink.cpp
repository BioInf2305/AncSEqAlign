#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <random>

/*
 *the following function collect the alternate bases with the quality score exceed by that of the defined by the user
 */

int countChunk = 0 ;

std::vector<std::string> collectAltBaseVec(std::vector<std::string>& vectorLine, int& qualScore){
	std::vector<std::string> qualAltAlleles;
	const char *altAlleles = vectorLine[3].c_str();
	const char *baseQual = vectorLine[4].c_str();
	for(int i=0;i<strlen(baseQual);i++){
	 	if(int(baseQual[i])>= qualScore+33){
			std::string highQualBase(2,altAlleles[i]);
			if(qualAltAlleles.size()==0){
				qualAltAlleles.push_back(highQualBase);
			}
			else {
			if(std::find(qualAltAlleles.begin(), qualAltAlleles.end(), highQualBase) == qualAltAlleles.end()){
				qualAltAlleles.push_back(highQualBase);
			}
		   }
	 	}
	}
	return qualAltAlleles;
}

/*
 *the following function write the vector to the ped and map file
 *
 */

void writeVectorToPlinkFiles(std::vector<std::vector<std::string>>& genoRecordVec, std::string& sampleName){
        std::string mapSu(".map");
	std::string pedSu(".ped");
        pedSu = sampleName+pedSu;
	mapSu = sampleName+mapSu;
        std::ofstream dest;
	dest.open(pedSu, std::ofstream::out | std::ofstream::app);
	std::ofstream dest1;
	dest1.open(mapSu, std::ofstream::out | std::ofstream::app);
	if (countChunk == 0){
	std::string famName(" 0 0 0 0");
	famName = sampleName+" "+sampleName+famName;
	dest<<famName;
	}
        for(int i=0; i<genoRecordVec.size();i++){
                        dest<<" "<<genoRecordVec[i][2][0]<<" "<<genoRecordVec[i][2][1];
			dest1<<genoRecordVec[i][0]<<" "<<genoRecordVec[i][0]<<"_"<<genoRecordVec[i][1]<<" 0 "<<genoRecordVec[i][1]<<"\n";
	}
	countChunk++;
} 

/*
 *
 *the following function read the file lines into vector
 */

void readFileToVector(const std::string& filename, int& qualScore, std::string& sampleName){
    std::ifstream source;

    std::string pedSu(".ped");

    //read pileup file
    
    source.open(filename);
    std::vector<std::vector<std::string>> lineVecs;
    std::string line;
    while (std::getline(source, line))
    {
	std::stringstream ss ( line );
	std::string word;
	std::vector<std::string> vectorLine;
	std::vector<std::string> outputLine;
	while ( std::getline ( ss, word, ' ' ) ){
		vectorLine.push_back(word);
	}
	if(vectorLine.size()!=5){std::cout<<"INPUT ERROR: the total numbers of column is not equal to five "<<line<<"\n";std::exit(EXIT_FAILURE);}
	std::vector<std::string> qualAltAlleles = collectAltBaseVec(vectorLine, qualScore);
	if(qualAltAlleles.size()<3 and qualAltAlleles.size()>0){
		outputLine.push_back(vectorLine[0]);
		outputLine.push_back(vectorLine[1]);
		if(qualAltAlleles.size()==2){
			srand((unsigned int)time(NULL));
			int num = rand()%2;
			outputLine.push_back(qualAltAlleles[num]);
			}
		else{
			outputLine.push_back(qualAltAlleles[0]);
		}
	lineVecs.push_back(outputLine);
		}
	//std::cout<<lineVecs.size()<<"\n";
	if(lineVecs.size() > 1000000){
		writeVectorToPlinkFiles(lineVecs, sampleName);
		std::cout<<" read another 1000000 positions "<<"\n";
		lineVecs.clear();
		}
	}
        std::ofstream dest;
	dest.open(pedSu, std::ofstream::out | std::ofstream::app);
	dest<<"\n";
    //return lineVecs;
}


int main(int argc, char **argv){
    std::string charactersFilename(argv[1]);
    int qualScore = std::stoi(argv[2]);
    std::string sampleName(argv[3]);
    readFileToVector(charactersFilename, qualScore, sampleName);
    //std::vector<std::vector<std::string>> genoRecordVec = readFileToVector(charactersFilename, qualScore);
    //writeVectorToPlinkFiles(genoRecordVec, sampleName);
}
