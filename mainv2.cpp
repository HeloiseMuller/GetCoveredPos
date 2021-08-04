/*
 * Copyright 2021, Héloïse Muller
 *
 * heloise.muller@egce.cnrs-gif.fr
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */

#include <vector>
#include <string>
#include <array>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdexcept>

#include <unistd.h>
#include <omp.h>


void outputIf(std::string chrom, int lpos, int memorize, int minLength, std::ofstream &file){
  if (memorize >= minLength){
    //if minLength ok, save this window in a bed format
      file << chrom << "\t" << 
      lpos  - memorize  << "\t" << //-1 because bed exclude first position
      lpos  << std::endl; //-1 because this a the pos < minDepth. We want the one before that
  }
} 

const char* OPTIONS = "h:c:d:l:o:";

void usage() {
  std::cout << " NameSofware -c <coverage file> -d <minimum depth> -l <minimum length> -o <output file> [-t <num threads>]  [-h]\n\n";

  std::cout << "   -c <coverage file>       : name of the file containing sequencing depth at each positions\n";
  std::cout << "   -d <minimum depthofile>  : minimum value to be considered as coverd. Default number is 1\n";
  std::cout << "   -l <minimum length>      : minimum length coverds, Default number is 300\n";
  std::cout << "   -o <output file>         : name of the file to save results in\n";
  std::cout << "   -h                       : Show this message and exit\n\n";

  std::cout << " Writen by Héloïse Muller (heloise.muller@egce.cnrs-gif.fr)\n";
  std::cout << " Copyright (c) 2021, Héloïse Muller\n";
  std::cout << " All rights reserved.\n\n";

  std::cout << " Released under the terms and conditions of the CeCILL-v2.1 license.\n";
}

int main(int argc, char** argv) {

  // Parameters
  std::string coverage_file;
  std::string out_file;
  int minDepth = 1;
  int minLength = 300;
  int nthreads = 1;

  // Parse command line arguments
  int opt = 0;
  while((opt = getopt(argc, argv, OPTIONS)) != -1) {
    switch (opt) {
      case 't':
        nthreads = std::stoi(optarg);
        break;

      case 'c':
        coverage_file = std::string(optarg);
        break;

      case 'd':
        minDepth = std::stoi(optarg);
        break;

      case 'l':
        minLength = std::stoi(optarg);
        break;

      case 'o':
        out_file = std::string(optarg);
        break;

      case 'h':
        usage();
        return 0;
        break;
         
      default:
        usage();
        return 1;
        break;
    }
  }

  if(coverage_file.size() == 0 ) {
    std::cerr << " ERROR: Improper arguments provided.\n\n";
    usage();
    return 1;
  }

  
  //Open stream to write in file
  std::ofstream ofile(out_file);
  std::ifstream ifile(coverage_file);

  //Variables for the three colomns
  std::string schrom, chrom;
  uint64_t pos;
  uint64_t lpos = 1; //Set last pos to 1
  uint64_t cov;
  int memorize = 0;
  int i = 0;

  //The name of the first chromosome is:
  ifile >> schrom;
  //Go back begining file
  ifile.clear();
  ifile.seekg(0);
 
  //Go through file but don't need to save it
  std::string element;
  while(ifile >> element) { //Go through file word by word
    if(i == 0) {
      chrom = element;
      i++;
    } else if(i == 1) {
      if(std::stoll(element) < 0) {
        std::cerr << "ERROR: position negative\n";
        std::exit(1);
      }
      pos = std::stoull(element); 
      i++;
    } else if(i == 2) {
      cov = std::stoull(element);
    
      // Reset i
      i = 0;

      //Check that we are still in the same chrom
      if (chrom == schrom){
        //count the position only if sequencing depth >= minDepth
        if (cov>=minDepth){
          memorize ++;
        } else {
          //when arrive at a position for which sequencing depth < minDepth, save the previous window if large enough
          outputIf(chrom, lpos, memorize, minLength, ofile);
          // reinitialize memorize
          memorize = 0;
        }
      }
    

      //We are on another chrom, so in a different window
      else {
        //Check if previous window large enough to be saved
        outputIf(schrom,  lpos, memorize, minLength, ofile);

        //Reinisialize for the new window
        schrom = chrom;
        if (cov>=minDepth){
          memorize = 1;
        } else {
          memorize = 0;
        } 
      }

      //Keep in mind the pos we just worked with (lpos)
      lpos = pos;

    }
  }
  
  //We read all the file. Check if the last window should be saved
  outputIf(chrom, lpos, memorize, minLength, ofile);
  ofile.close();
   
    if (i != 0 && element != "\n") {
    std::cerr << " ERROR: Unexpected dangling values in windows file\n";
    std::exit(1);
  }

  ifile.close();


  return 0;
}
