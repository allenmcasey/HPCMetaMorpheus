#include "CrosslinkSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "Crosslinker.h"

#include "Sort.h"
#include "BinaryPack.h"

#include <numeric>
#include <algorithm>
#include <sstream>
#include <iomanip>

using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
    namespace CrosslinkSearch
    {
        
        CrosslinkSpectralMatch::CrosslinkSpectralMatch(PeptideWithSetModifications *theBestPeptide, int notch,
                                                       double score, int scanIndex,
                                                       Ms2ScanWithSpecificMass *scan,
                                                       DigestionParams *digestionParams,
                                                       std::vector<MatchedFragmentIon*> &matchedFragmentIons) :
            PeptideSpectralMatch(theBestPeptide, notch, score, scanIndex, scan, digestionParams, matchedFragmentIons)
        {
            this->setXLTotalScore(score);
        }

        CrosslinkSpectralMatch::CrosslinkSpectralMatch(PeptideWithSetModifications *theBestPeptide, int notch,
                                                       double score, int scanIndex,
                                                       std::string scanFullFilePath, int scanOneBasedScanNumber,
                                                       std::optional<int> scanOneBasedPrecursorScanNumber,
                                                       double scanRetentionTime, int scanNumPeaks, double scanTotalIonCurrent,
                                                       int scanPrecursorCharge, double scanPrecursorMonoisotopicPeakMz,
                                                       double scanPrecursorMass,                                                       
                                                       DigestionParams *digestionParams,
                                                       std::vector<MatchedFragmentIon*> &matchedFragmentIons) :
            PeptideSpectralMatch(theBestPeptide, notch, score, scanIndex, scanFullFilePath,
                                 scanOneBasedScanNumber, scanOneBasedPrecursorScanNumber, scanRetentionTime,
                                 scanNumPeaks, scanTotalIonCurrent, scanPrecursorCharge, scanPrecursorMonoisotopicPeakMz,
                                 scanPrecursorMass, digestionParams, matchedFragmentIons)
        {
            this->setXLTotalScore(score);
        }

        
        CrosslinkSpectralMatch *CrosslinkSpectralMatch::getBetaPeptide() const
        {
            return privateBetaPeptide;
        }
        
        void CrosslinkSpectralMatch::setBetaPeptide(CrosslinkSpectralMatch *value)
        {
            privateBetaPeptide = value;
        }
        
        std::vector<int> CrosslinkSpectralMatch::getLinkPositions() const
        {
            return privateLinkPositions;
        }
        
        void CrosslinkSpectralMatch::setLinkPositions(const std::vector<int> &value)
        {
            privateLinkPositions = value;
        }
        
        double CrosslinkSpectralMatch::getDeltaScore() const
        {
            return privateDeltaScore;
        }
        
        void CrosslinkSpectralMatch::setDeltaScore(double value)
        {
            privateDeltaScore = value;
        }
        
        double CrosslinkSpectralMatch::getXLTotalScore() const
        {
            return privateXLTotalScore;
        }
        
        void CrosslinkSpectralMatch::setXLTotalScore(double value)
        {
            privateXLTotalScore = value;
        }
        
        int CrosslinkSpectralMatch::getXlProteinPos() const
        {
            return privateXlProteinPos;
        }
        
        void CrosslinkSpectralMatch::setXlProteinPos(int value)
        {
            privateXlProteinPos = value;
        }
        
        std::vector<int> CrosslinkSpectralMatch::getXlRank() const
        {
            return privateXlRank;
        }
        
        void CrosslinkSpectralMatch::setXlRank(const std::vector<int> &value)
        {
            privateXlRank = value;
        }
        
        std::string CrosslinkSpectralMatch::getParentIonExist() const
        {
            return privateParentIonExist;
        }
        
        void CrosslinkSpectralMatch::setParentIonExist(const std::string &value)
        {
            privateParentIonExist = value;
        }
        
        int CrosslinkSpectralMatch::getParentIonExistNum() const
        {
            return privateParentIonExistNum;
        }
        
        void CrosslinkSpectralMatch::setParentIonExistNum(int value)
        {
            privateParentIonExistNum = value;
        }
        
        std::vector<int> CrosslinkSpectralMatch::getParentIonMaxIntensityRanks() const
        {
            return privateParentIonMaxIntensityRanks;
        }
        
        void CrosslinkSpectralMatch::setParentIonMaxIntensityRanks(const std::vector<int> &value)
        {
            privateParentIonMaxIntensityRanks = value;
        }
        
        PsmCrossType CrosslinkSpectralMatch::getCrossType() const
        {
            return privateCrossType;
        }
        
        void CrosslinkSpectralMatch::setCrossType(PsmCrossType value)
        {
            privateCrossType = value;
        }
        
        std::vector<int> CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(std::vector<char> &crosslinkerModSites,
                                                                                PeptideWithSetModifications *peptide)
        {
            std::vector<int> possibleXlPositions;

            bool wildcard =false;
            for ( char p: crosslinkerModSites ){
                if (  p == 'X' ) {
                    wildcard = true;
                    break;
                }                      
            }
            
            for (int r = 0; r < (int) peptide->getBaseSequence().size(); r++)
            {
                //if (crosslinkerModSites.Contains(peptide->getBaseSequence()[r]) || wildcard)
                if ( std::find(crosslinkerModSites.begin(), crosslinkerModSites.end(), peptide->getBaseSequence()[r]) !=
                     crosslinkerModSites.end() || wildcard ) 
                {
                    possibleXlPositions.push_back(r + 1);
                }
            }
            
            return possibleXlPositions;
        }
        
        std::vector<int> CrosslinkSpectralMatch::GenerateIntensityRanks(std::vector<double> &experimental_intensities)
        {
            auto y = experimental_intensities;
#ifdef ORIG
            auto x = Enumerable::Range(1, y.size()).OrderBy([&] (std::any p) {
                    return p;
                })->ToArray();
            Array::Sort(y, x);
#endif
            std::vector<int> x(y.size());
            std::iota(x.begin(), x.end(), 1);
            Sort::SortPairs(y, x, y.size() );
            
#ifdef ORIG
            auto experimental_intensities_rank = Enumerable::Range(1, y.size()).OrderByDescending([&] (std::any p)  {
                    return p;
                })->ToArray();
            Array::Sort(x, experimental_intensities_rank);
#endif
            std::vector<int> experimental_intensities_rank(y.size() );
            int n = y.size();
            std::generate(experimental_intensities_rank.begin(), experimental_intensities_rank.end(), [&] () {return n--;});
            Sort::SortPairs(x, experimental_intensities_rank, x.size() );
            
            return experimental_intensities_rank;
        }
        
        std::string CrosslinkSpectralMatch::GetTabSepHeaderCross()
        {
            auto sb = new StringBuilder();
            sb->append("File Name\t");
            sb->append("Scan Number\t");
            sb->append("Precursor Scan Number\t");
            sb->append("Precursor MZ\t");
            sb->append("Precursor Charge\t");
            sb->append("Precursor Mass\t");
            sb->append("Cross Type\t");
            sb->append("Link Residues\t");
            
            sb->append("Peptide\t");
            sb->append("Protein Accession\t");
            sb->append("Protein Link Site\t");
            sb->append("Base Sequence\t");
            sb->append("Full Sequence\t");
            sb->append("Peptide Monoisotopic Mass\t");
            sb->append("Score\t");
            sb->append("Rank\t");
            
            sb->append("Matched Ion Series\t");
            sb->append("Matched Ion Mass-To-Charge Ratios\t");
            sb->append("Matched Ion Mass Diff (Da)\t");
            sb->append("Matched Ion Mass Diff (Ppm)\t");
            sb->append("Matched Ion Intensities\t");
            sb->append("Matched Ion Counts\t");
            
            sb->append("Beta Peptide\t");
            sb->append("Beta Peptide Protein Accession\t");
            sb->append("Beta Peptide Protein LinkSite\t");
            sb->append("Beta Peptide Base Sequence\t");
            sb->append("Beta Peptide Full Sequence\t");
            sb->append("Beta Peptide Theoretical Mass\t");
            sb->append("Beta Peptide Score\t");
            sb->append("Beta Peptide Rank\t");
            
            sb->append("Beta Peptide Matched Ion Series\t");
            sb->append("Beta Peptide Matched Ion Mass-To-Charge Ratios\t");
            sb->append("Beta Peptide Matched Ion Mass Diff (Da)\t");
            sb->append("Beta Peptide Matched Ion Mass Diff (Ppm)\t");
            sb->append("Beta Peptide Matched Ion Intensities\t");
            sb->append("Beta Peptide Matched Ion Counts\t");
            
            sb->append("Summary\t");
            sb->append("XL Total Score\t");
            sb->append("Mass Diff (Da)\t");
            sb->append("Parent Ions\t");
            sb->append("ParentIonsNum\t");
            sb->append("ParentIonMaxIntensityRank\t");
            sb->append("Decoy/Contaminant/Target\t");
            sb->append("QValue\t");
            
            
            std::string s = sb->toString();
            delete sb;
            return s;
        }
        
        std::string CrosslinkSpectralMatch::GetTabSepHeaderSingle()
        {
            auto sb = new StringBuilder();
            sb->append("File Name\t");
            sb->append("Scan Number\t");
            sb->append("Precursor Scan Number\t");
            sb->append("Precursor MZ\t");
            sb->append("Precursor Charge\t");
            sb->append("Precursor Mass\t");
            sb->append("Cross Type\t");
            sb->append("Link Residues\t");

            sb->append("Protein Accession\t");
            sb->append("Protein Link Site\t");
            sb->append("Base Sequence\t");
            sb->append("Full Sequence\t");
            sb->append("Peptide Monoisotopic Mass\t");
            sb->append("Score\t");
            sb->append("Rank\t");
            
            sb->append("Matched Ion Series\t");
            sb->append("Matched Ion Mass-To-Charge Ratios\t");
            sb->append("Matched Ion Mass Diff (Da)\t");
            sb->append("Matched Ion Mass Diff (Ppm)\t");
            sb->append("Matched Ion Intensities\t");
            sb->append("Matched Ion Counts\t");
            sb->append("Decoy/Contaminant/Target\t");
            sb->append("QValue\t");
            
            std::string s = sb->toString();
            delete sb;
            return s;
        }
        
        std::string CrosslinkSpectralMatch::GetTabSepHeaderGlyco()
        {
            auto sb = new StringBuilder();
            sb->append("File Name\t");
            sb->append("Scan Number\t");
            sb->append("Precursor Scan Number\t");
            sb->append("Precursor MZ\t");
            sb->append("Precursor Charge\t");
            sb->append("Precursor Mass\t");
            sb->append("Cross Type\t");
            sb->append("Link Residues\t");

            sb->append("Protein Accession\t");
            sb->append("Protein Link Site\t");
            sb->append("Base Sequence\t");
            sb->append("Full Sequence\t");
            sb->append("Peptide Monoisotopic Mass\t");
            sb->append("Score\t");
            sb->append("Rank\t");
            
            sb->append("Matched Ion Series\t");
            sb->append("Matched Ion Mass-To-Charge Ratios\t");
            sb->append("Matched Ion Mass Diff (Da)\t");
            sb->append("Matched Ion Mass Diff (Ppm)\t");
            sb->append("Matched Ion Intensities\t");
            sb->append("Matched Ion Counts\t");
            
            sb->append("Decoy/Contaminant/Target\t");
            sb->append("QValue\t");
            
            sb->append("GlyID\t");
            sb->append("GlyMass\t");
            sb->append("GlyStruct(H,N,A,G,F)\t");
            
            std::string s = sb->toString();
            delete sb;
            return s;
        }
        
        std::string CrosslinkSpectralMatch::ToString()
        {
            std::string position = "";
            switch (getCrossType())
            {
                case PsmCrossType::Single:
                    break;
                    
                case PsmCrossType::Loop:
                    position = "(" + std::to_string(getLinkPositions()[0]) + "-" + std::to_string(getLinkPositions()[1]) + ")";
                    break;
                    
                default:
                    position = "(" + std::to_string(getLinkPositions()[0]) + ")";
                    break;
            }
            
            auto sb = new StringBuilder();
            sb->append(getFullFilePath() + "\t");
            sb->append(std::to_string(getScanNumber()) + "\t");
            if ( getPrecursorScanNumber().has_value() ) {
                sb->append(std::to_string(getPrecursorScanNumber().value()) + "\t");
            }
            else {
                std::string s = "-\t";
                sb->append(s);
            }
            sb->append(std::to_string(getScanPrecursorMonoisotopicPeakMz()) + "\t");
            sb->append(std::to_string(getScanPrecursorCharge()) + "\t");
            sb->append(std::to_string(getScanPrecursorMass()) + "\t");
            auto crosslinktype = getCrossType();
            sb->append(PsmCrossTypeToString(crosslinktype) + "\t");
            
            if (getLinkPositions().size() > 0)
            {
                if (getCrossType() == PsmCrossType::Loop)
                {
                    std::stringstream ss;
                    ss << getBaseSequence()[getLinkPositions()[0] - 1] << ';' <<
                        getBaseSequence()[getLinkPositions()[1] - 1] << '\t';
                    sb->append(ss.str() );
                }
                else if (getCrossType() == PsmCrossType::Inter || getCrossType() == PsmCrossType::Intra)
                {
                    std::stringstream ss;
                    ss << getBaseSequence()[getLinkPositions()[0] - 1] << ';' <<
                        getBetaPeptide()->getBaseSequence()[getBetaPeptide()->getLinkPositions()[0] - 1] << '\t';
                    sb->append(ss.str());
                }
                else
                {
                    // deadend
                    std::stringstream ss;
                    ss << getBaseSequence()[getLinkPositions()[0] - 1] << "\t";
                    sb->append(ss.str() );
                }
            }
            else
            {
                sb->append("\t");
            }
            
            sb->append("\t");
            sb->append(getProteinAccession() + "\t");
            sb->append(std::to_string(getXlProteinPos()) + "\t");
            sb->append(getBaseSequence() + "\t");
            sb->append(getFullSequence() + position + "\t");
            sb->append((getPeptideMonisotopicMass().has_value() ? std::to_string(getPeptideMonisotopicMass().value()) : "---"));
            sb->append("\t");
            sb->append(std::to_string(getScore()) + "\t");
            sb->append(std::to_string(getXlRank()[0]) + "\t");
            
            for (auto mid : MatchedIonDataDictionary(this))
            {
                sb->append(std::get<1>(mid));
                sb->append("\t");
            }
            
            if (getBetaPeptide() != nullptr)
            {
                auto betaPeptide = getBetaPeptide();
                
                sb->append("\t");
                sb->append(betaPeptide->getProteinAccession() + "\t");
                sb->append(std::to_string(betaPeptide->getXlProteinPos()) + "\t");
                sb->append(betaPeptide->getBaseSequence() + "\t");
                sb->append(betaPeptide->getFullSequence() + "(" + std::to_string(betaPeptide->getLinkPositions()[0]) +
                           ")" + "\t");
                sb->append(std::to_string(betaPeptide->getPeptideMonisotopicMass().value()) + "\t");
                sb->append(std::to_string(betaPeptide->getScore()) + "\t");
                sb->append(std::to_string(getXlRank()[1]) + "\t");
                
                for (auto betamid : MatchedIonDataDictionary(this->getBetaPeptide()))
                {
                    sb->append(std::get<1>(betamid));
                    sb->append("\t");
                }
                
                sb->append("\t");
                sb->append(std::to_string(getXLTotalScore()) + "\t");
                
                // mass of crosslinker
                sb->append(((getPeptideMonisotopicMass().has_value()) ? std::to_string(getScanPrecursorMass() -
                            betaPeptide->getPeptideMonisotopicMass().value() - getPeptideMonisotopicMass().value()) : "---"));
                sb->append("\t");
                
                int alphaNumParentIons = 0;
                for ( auto p : getMatchedFragmentIons() ) {
                    if ( p->NeutralTheoreticalProduct->productType == ProductType::M ) {
                        alphaNumParentIons++;
                    }
                }

                int betaNumParentIons = 0;
                for ( auto p :  betaPeptide->getMatchedFragmentIons() ) {
                    if ( p->NeutralTheoreticalProduct->productType == ProductType::M ) {
                        alphaNumParentIons++;
                    }
                }

                
                sb->append(std::to_string(alphaNumParentIons) + ";" + std::to_string(betaNumParentIons) + "\t");
                sb->append(std::to_string(alphaNumParentIons) + std::to_string(betaNumParentIons) + "\t");
                sb->append(((getParentIonMaxIntensityRanks().size() > 0) && (!getParentIonMaxIntensityRanks().empty()) ?
                   std::to_string(*std::min_element(getParentIonMaxIntensityRanks().begin(), getParentIonMaxIntensityRanks().end()))
                            : "-"));
                sb->append("\t");                            
            }
            
            if (getBetaPeptide() == nullptr)
            {
                sb->append((getIsDecoy()) ? "D" : (getIsContaminant()) ? "C" : "T");
                sb->append("\t");
            }
            else
            {
                sb->append((getIsDecoy() || getBetaPeptide()->getIsDecoy()) ? "D" :
                           (getIsContaminant() || getBetaPeptide()->getIsContaminant()) ? "C" : "T");
                sb->append("\t");
            }
            
            sb->append(std::to_string(getFdrInfo()->getQValue()));
            sb->append("\t");
            
            
            std::string s= sb->toString();
            delete sb;
            return s;
        }
        
        std::vector<std::tuple<std::string, std::string>> CrosslinkSpectralMatch::MatchedIonDataDictionary(PeptideSpectralMatch *psm)
        {
            std::vector<std::tuple<std::string, std::string>> s;
            AddMatchedIonsData(s, psm);
            return s;
        }

        /**==================================================/
         *                   NEW VERSION                    *
         * =================================================*/

        int CrosslinkSpectralMatch::Pack(char *buf, size_t &buf_len, CrosslinkSpectralMatch *csm)
        {
            std::vector<CrosslinkSpectralMatch *> csmVec;
            csmVec.push_back(csm);
            int pos = CrosslinkSpectralMatch::Pack(buf, buf_len, csmVec);
            
            return pos;            
        }

        int CrosslinkSpectralMatch::Pack(char *buf, size_t &buf_len, const std::vector<CrosslinkSpectralMatch *> &csmVec)
        {
	    std::cout << "Packing" << std::endl;
            std::vector<SerializedCrosslinkSpectralMatch> serializedCSMVec;

            for (auto csm : csmVec)
            {
                SerializedCrosslinkSpectralMatch serializedCSM = CrosslinkSpectralMatch::Pack_internal(csm);
                serializedCSMVec.push_back(serializedCSM);

		std::cout << "\tpushed to vector..." << std::endl;

                auto betaPeptide = csm->getBetaPeptide();

		std::cout << "\tChecked for betapep..."  << std::endl;		

                if (betaPeptide != nullptr)
                {
		    std::cout << "\t\thas betapep..."  << std::endl;
                    SerializedCrosslinkSpectralMatch serializedBetaPeptide = CrosslinkSpectralMatch::Pack_internal(betaPeptide);
                    serializedCSMVec.push_back(serializedBetaPeptide);
                }
		std::cout << "\tCSM complete..."  << std::endl;
            }

	    std::cout << "Finished internal packing, vector size: " << serializedCSMVec.size() << std::endl;

            std::stringstream sstream;
            msgpack::pack(sstream, serializedCSMVec);
            std::string tmp = sstream.str();            
            int bufSize = sizeof(tmp);

	    std::cout << "Size of buffer: " << bufSize << std::endl;

	    memcpy(buf, tmp.c_str(), bufSize);
	    std::cout << buf << std::endl;
            return bufSize;
        }

        SerializedCrosslinkSpectralMatch CrosslinkSpectralMatch::Pack_internal(CrosslinkSpectralMatch *csm)
        {
	    std::cout << "Packing internally..." << std::endl;
            std::vector<int> lPositions  = csm->getLinkPositions();
            std::vector<int> xlRanks = csm->getXlRank();

	    std::cout << "1"  << std::endl;

            bool hasNotchValue = csm->getNotch().has_value();
            int notchValue = 0;
            if (hasNotchValue)
            {
                notchValue = csm->getNotch().value();
            }

	    std::cout << "2"  << std::endl;

            double xlTotalScore = csm->getXLTotalScore();
            double deltaScore = csm->getDeltaScore();
            double score = csm->getScore();
            double runnerUpScore = csm->getRunnerUpScore();
            double peptideMonoisotopicMass = csm->getPeptideMonisotopicMass().value();

	    std::cout << "3"  << std::endl;

            int scanNumber = csm->getScanNumber();
            int xlProteinPos = csm->getXlProteinPos();
            int matchedFragmentIonsSize = (int)csm->getMatchedFragmentIons().size();
            int lPositionsSize = (int)lPositions.size();
            int xlRanksSize = (int)xlRanks.size();

            bool hasBetaPeptide = csm->getBetaPeptide() != nullptr; 

	    std::cout << "4"  << std::endl;

            PsmCrossType ctype = csm->getCrossType();
            std::string psmCrossTypeString = PsmCrossTypeToString(ctype);
            
	    std::cout << "5"  << std::endl;

            auto tvar = csm->getPrecursorScanNumber();
            bool hasPrecursorScanNumber = tvar.has_value();
            int precursorScanNumber = 0;
            if (hasPrecursorScanNumber)
            {
                precursorScanNumber = tvar.value();
            }

	    std::cout << "6"  << std::endl;

            int scanExperimentalPeaks = csm->getScanExperimentalPeaks();
            int scanPrecursorCharge = csm->getScanPrecursorCharge();

	    std::cout << "7"  << std::endl;

            double scanRetentionTime = csm->getScanRetentionTime();
            double totalIonCurrent = csm->getScanRetentionTime();
            double scanPrecursorMonoisotopicPeakMz = csm->getScanPrecursorMonoisotopicPeakMz();
            double scanPrecursorMass = csm->getScanPrecursorMass();

            std::string fullFilePath = csm->getFullFilePath();

	    std::cout << "8"  << std::endl;

            FdrInfo *fdr = csm->getFdrInfo();
	    std::cout << "Got fdr new"  << std::endl;

            bool has_fdr = (fdr != nullptr);

	    std::cout << has_fdr  << std::endl;

	    std::cout << "sumpin" << std::endl;

	    SerializedFdrInfo sFdr;	
	    if (has_fdr)
            {

		std::cout << "serializaing fdr..." << std::endl;
	
                sFdr = {
                    fdr->getCumulativeTarget(),
                    fdr->getCumulativeDecoy(),
                    fdr->getQValue(),
                    fdr->getCumulativeTargetNotch(),
                    fdr->getCumulativeDecoyNotch(),
                    fdr->getQValueNotch(),
                    fdr->getMaximumLikelihood(),
                    fdr->getEValue(),
                    fdr->getEScore(),
                    fdr->getCalculateEValue(),
                    has_fdr,
                };
            }
            else
            {
		std::cout << "not serializaing fdr because it doesn't exist..." << std::endl;

                sFdr = {
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                    NULL, NULL, false,
                };
            }		

 	    std::cout << "Done serializing fdr"  << std::endl;

            auto dp = csm->digestionParams;
            std::string digestionParamsString = dp->ToString();

            auto uMapPep = csm->getPeptidesToMatchingFragments(); 
            if ( uMapPep.size() != 1 ) {
                std::cout << "CrosslinkSpectralMatch::Pack: Error - unordered_map has more than one entry!\n";
            }
            auto pep = std::get<0>(*uMapPep.begin());

            SerializedPeptide sPep = {
                pep->getOneBasedStartResidueInProtein(),
                pep->getOneBasedEndResidueInProtein(),
                pep->getMissedCleavages(),
                pep->NumFixedMods,
                pep->getPeptideDescription(),
                CleavageSpecificityExtension::GetCleavageSpecificityAsString(pep->getCleavageSpecificityForFdrCategory()),
                pep->getFullSequence(),
                pep->getDigestionParamString(),
            };

	    std::cout << "Done serilizaing SPep"  << std::endl;

            std::string accession = pep->getProteinAccession();
            if ( accession != "" )  {
                sPep.GetProteinAccession = accession;                                     
            }
            else  {
                sPep.GetProteinAccession = "-";
            }

            std::vector<SerializedMatchedFragmentIon> matchedFragmentIons;
            auto mFrIons = csm->getMatchedFragmentIons();
            for (auto ion : mFrIons) 
            {
		auto product = ion->NeutralTheoreticalProduct;
                auto tfragment = product->TerminusFragment;
	        auto ptype = product->productType;

                SerializedMatchedFragmentIon sIon = {
                    ion->Mz,
                    ion->Intensity,
                    ion->Charge,
                    product->NeutralLoss,
                    Fragmentation::ProductTypeToString(ptype),
                    Fragmentation::FragmentationTerminusToString(tfragment->Terminus),
                    tfragment->NeutralMass,
                    tfragment->FragmentNumber,
                    tfragment->AminoAcidPosition,
                };

                matchedFragmentIons.push_back(sIon);
            }

	    std::cout << "Done serializing MFR"  << std::endl;

            SerializedCrosslinkSpectralMatch serializedCSM = {
                hasNotchValue, 
                notchValue, 
                xlTotalScore, 
                deltaScore, 
                score, 
                runnerUpScore,
                peptideMonoisotopicMass, 
                scanNumber, 
                xlProteinPos, 
                matchedFragmentIonsSize,
                lPositionsSize, 
                xlRanksSize, 
                hasBetaPeptide, 
                psmCrossTypeString, 
                hasPrecursorScanNumber,
                precursorScanNumber, 
                scanExperimentalPeaks, 
                scanPrecursorCharge, 
                scanRetentionTime,
                totalIonCurrent, 
                scanPrecursorMonoisotopicPeakMz, 
                scanPrecursorMass, 
                fullFilePath,
                sFdr, 
                lPositions, 
                xlRanks, 
                digestionParamsString, 
                sPep, 
                matchedFragmentIons,
            };

	    std::cout << "Done internally"  << std::endl;

            return serializedCSM;
        }

        void CrosslinkSpectralMatch::Unpack (char *buf, size_t buf_len, int count, size_t &len,
                                             std::vector<CrosslinkSpectralMatch*> &pepVec,
                                             const std::vector<Modification*> &mods,
                                             const std::vector<Protein *> &proteinList)
        {
            std::stringstream sstream;
            sstream << buf;
	    std::cout << "sstream size before unpacking: " << sstream.str().size() << std::endl;

            msgpack::object_handle objHandle = msgpack::unpack(sstream.str().data(), sstream.str().size());
            msgpack::object const& obj = objHandle.get();
            auto serializedCSMVec = obj.as< std::vector<SerializedCrosslinkSpectralMatch> >();

            int index = 0;
            while (index < serializedCSMVec.size())
            {
                CrosslinkSpectralMatch *pep;
                bool has_beta_peptide = false;

                CrosslinkSpectralMatch::Unpack_internal (serializedCSMVec[index++], &pep, mods, proteinList, has_beta_peptide);
                pepVec.push_back(pep);

                if ( has_beta_peptide ) {
                    CrosslinkSpectralMatch *beta_pep;
		    has_beta_peptide = false;
                    CrosslinkSpectralMatch::Unpack_internal(serializedCSMVec[index++], &beta_pep, mods, proteinList, has_beta_peptide);
                    pep->setBetaPeptide(beta_pep);
                }
            }
        }

        void CrosslinkSpectralMatch::Unpack (char *buf, size_t buf_len, size_t &len,
                                             CrosslinkSpectralMatch** newCsm,
                                             const std::vector<Modification*> &mods,
                                             const std::vector<Protein *> &proteinList )
        {
            std::stringstream sstream;
            sstream << buf;

            msgpack::object_handle objHandle = msgpack::unpack(sstream.str().data(), sstream.str().size());
            msgpack::object const& obj = objHandle.get();
            auto serializedCSMVec = obj.as< std::vector<SerializedCrosslinkSpectralMatch> >();

            bool has_beta_peptide=false;            
            CrosslinkSpectralMatch::Unpack_internal (serializedCSMVec[0], newCsm, mods, proteinList, has_beta_peptide);
            if ( has_beta_peptide) 
            {
                CrosslinkSpectralMatch* beta_pep;

                CrosslinkSpectralMatch::Unpack_internal (serializedCSMVec[1], &beta_pep, mods, proteinList, has_beta_peptide);
                (*newCsm)->setBetaPeptide(beta_pep);
            }            
        }

        void CrosslinkSpectralMatch::Unpack_internal (SerializedCrosslinkSpectralMatch sCSM, CrosslinkSpectralMatch **newCSM,
                                                    const std::vector<Modification*> &mods, const std::vector<Protein *> &proteinList,
                                                    bool &has_beta_peptide )
        { 
	    int notch = -1, scannumber, proteinPos, matchedFragmentIonsVecsize, lpositionsize, xlranksize;
            double  deltaScore, XLTotalScore, score, runnerUpScore, peptideMonisotopicMass = 0;
            bool  tmpvar;
            
            if ( sCSM.hasNotchValue) {
                notch = sCSM.notchValue;
            }

            XLTotalScore = sCSM.xlTotalScore;
            deltaScore = sCSM.deltaScore;
            score = sCSM.score;
            runnerUpScore = sCSM.runnerUpScore;
            peptideMonisotopicMass = sCSM.peptideMonoisotopicMass;

            scannumber = sCSM.scanNumber;
            proteinPos = sCSM.xlProteinPos;
            matchedFragmentIonsVecsize = sCSM.matchedFragmentIonsSize;
            lpositionsize = sCSM.lPositionsSize;
            xlranksize = sCSM.xlRanksSize;

            has_beta_peptide = sCSM.hasBetaPeptide;

            PsmCrossType ctype = PsmCrossTypeFromString(sCSM.psmCrossTypeToString);
            
            bool has_tvar = sCSM.hasPrecursorScanNumber;
            int scanPrecursorScanNumber = -1;
            if (has_tvar)
            {
                scanPrecursorScanNumber = sCSM.precursorScanNumber;
            }

            int scanExperimentalPeaks = sCSM.scanExperimentalPeaks;
            int scanPrecursorCharge = sCSM.scanPrecursorCharge;

            double scanRetentionTime = sCSM.scanRetentionTime;
            double scanTotalIonCurrent = sCSM.totalIonCurrent;
            double scanPrecursorMonoisotopicPeakMz = sCSM.scanPrecursorMonoisotopicPeakMz;
            double scanPrecursorMass = sCSM.scanPrecursorMass;

            std::string scanFullFilePath = sCSM.fullFilePath;
            
            //line 2: FdrInfo related data
            FdrInfo* fdr = nullptr;
            FdrInfo* tempFdr = new FdrInfo();
		
	    SerializedFdrInfo sFdrInfo = sCSM.fdrInfo;
            tempFdr->setCumulativeTarget(sFdrInfo.cumulativeTarget);
            tempFdr->setCumulativeDecoy(sFdrInfo.cumulativeDecoy);
            tempFdr->setQValue(sFdrInfo.qValue);
            tempFdr->setCumulativeTargetNotch(sFdrInfo.cumulativeTargetNotch);
            tempFdr->setCumulativeDecoyNotch(sFdrInfo.cumulativeDecoyNotch);
            tempFdr->setQValueNotch(sFdrInfo.qValueNotch);
            tempFdr->setMaximumLikelihood(sFdrInfo.maximumLikelihood);
            tempFdr->setEScore(sFdrInfo.eScore);
            tempFdr->setEValue(sFdrInfo.eValue);
            tempFdr->setCalculateEValue(sFdrInfo.calculateEValue);

            fdr = tempFdr;

            //line 3: linkPositions          
            std::vector<int> linkPosvec = sCSM.lPositions;
            
            //line 4: xlRank
            std::vector<int> xlRankVec = sCSM.xlRanks;
            
            //line 5: DigestionParams
            std::string dpstring = sCSM.digestionParamsString;
            DigestionParams *dp = DigestionParams::FromString(dpstring);

            /*******************************************************
             * Everything up to this point is good to go it seems...
             * proceed to checking out the PepWithSetMods & MFI
            */
		
	    SerializedPeptide sPep = sCSM.peptide;

            //line 6-10: PeptideWithSetModification
            CleavageSpecificity cvs = CleavageSpecificityExtension::ParseString(sPep.GetCleavageSpecificityAsString);

            std::unordered_map<std::string, Modification*> umsM;
            auto pep = new PeptideWithSetModifications(sPep.GetFullSequence, 
                                                          umsM,                
                                                          sPep.NumFixedMods,
                                                          dp,                  
                                                          nullptr,  
                                                          sPep.GetOneBasedStartResidueInProtein,
                                                          sPep.GetOneBasedEndResidueInProtein,
                                                          sPep.GetMissedCleavages,
                                                          cvs,
                                                          sPep.GetPeptideDescription);
            pep->SetNonSerializedPeptideInfo(mods, proteinList);	
	
            // protein accession?

#ifdef DEBUG
            // Safety check:
            double monIsotopicMass = pep->getMonoisotopicMass();
            if (monIsotopicMass != peptideMonisotopicMass) {
                std::cout << "Safety check failed when reconstructing PeptideWithSetMOdifications in CrosslinkSpectralMatch::Unpack(). " <<
                    "monIsotopicMass is " << monIsotopicMass << " should be " << peptideMonisotopicMass << std::endl;
            }
#endif
            
            // line 11-x: Vector of MatchedFragmentIons
            std::vector<MatchedFragmentIon*> matchedFragmentIonsVec;
            std::vector<SerializedMatchedFragmentIon> serializedMatchedFragmentIons = sCSM.matchedFragmentIons;

            for (auto i = 0; i < matchedFragmentIonsVecsize; i++) {

                SerializedMatchedFragmentIon sMaFObject = serializedMatchedFragmentIons[i];

                // get Ion params
                double mz = sMaFObject.mz;
                double intensity = sMaFObject.intensity;
                int charge = sMaFObject.charge;

                // get Product params
                double neutralLoss = sMaFObject.neutralLoss;
                auto pType = Fragmentation::ProductTypeFromString(sMaFObject.productType);

                // get NeutralTerminusFragment params
                double neutralMass = sMaFObject.neutralMass;
                auto fragmentationTerminus = Fragmentation::FragmentationTerminusFromString(sMaFObject.fragmentationTerminus);
                int fragmentNumber = sMaFObject.fragmentNumber;
                int aminoAcidPosition = sMaFObject.aminoAcidPosition;

                // create MatchedFragmentIon
                auto terminusFragment = new NeutralTerminusFragment(fragmentationTerminus, neutralMass, fragmentNumber, aminoAcidPosition);
                auto product = new Product(pType, terminusFragment, neutralLoss);
                auto newMaF = new MatchedFragmentIon(product, mz, intensity, charge);

                // add Ion to result vector
                matchedFragmentIonsVec.push_back(newMaF);
            }

            
            // We are trearint scannumber and scanindex as the same here. First, it is actually really the same
            // in many scenarios. Second, scanindex is not really used for the subsequent operations as far
            // as I can see right now.
            CrosslinkSpectralMatch *csm = new CrosslinkSpectralMatch (pep, notch, XLTotalScore, scannumber,
                                                                        scanFullFilePath, scannumber,
                                                                        scanPrecursorScanNumber,
                                                                        scanRetentionTime, scanExperimentalPeaks,
                                                                        scanTotalIonCurrent, scanPrecursorCharge,
                                                                        scanPrecursorMonoisotopicPeakMz,
                                                                        scanPrecursorMass, dp,
                                                                        matchedFragmentIonsVec);
            csm->setXLTotalScore(XLTotalScore);
            csm->setDeltaScore(deltaScore);
            csm->setXlProteinPos(proteinPos);
            csm->setScore(score);
            csm->setRunnerUpScore(runnerUpScore);
            csm->setCrossType (ctype);
            csm->setXlRank(xlRankVec);
            csm->setLinkPositions(linkPosvec);
            if (fdr != nullptr) 
            {
                csm->setFdrInfo(fdr);
            }
            csm->ResolveAllAmbiguities();

            // This is a hack. Need to find a proper solution.
            // In some situations involving Mods, the PeptideMonoisotopicMass
            // is not correct after deserialization.
            csm->setPeptideMonisotopicMass(peptideMonisotopicMass);

            *newCSM = csm;

            return ;
        }

    }
}
