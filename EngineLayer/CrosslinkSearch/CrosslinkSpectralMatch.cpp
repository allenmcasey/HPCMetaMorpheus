#include "CrosslinkSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "Crosslinker.h"

#include "Sort.h"
#include "BinaryPack.h"
#include "CrosslinkSpectralMatch.ph.h"

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

        // done
        int CrosslinkSpectralMatch::Pack(char *buf, size_t &buf_len, const std::vector<CrosslinkSpectralMatch *> &csmVec)
        {
            SerializedCrosslinkSpectralMatchVec sCsmVec;

            for (auto csm: csmVec) {
                
                // pack CSM and push to vector
                SerializedCrosslinkSpectralMatch sCsm = CrosslinkSpectralMatch::Pack_internal(csm);
                sCsmVec.push_back(sCsm);
                
                // check if csm has betapep
                auto betaPeptide = csm->getBetaPeptide();
                if (betaPeptide != nullptr) {

                    // pack betapep and push to vector
                    SerializedCrosslinkSpectralMatch sBetaPeptide = CrosslinkSpectralMatch::Pack_internal(betaPeptide);
                    sCsmVec.push_back(sBetaPeptide);                    
                }
            }

            // serialize csm vec
            std::string dataString;
            sCsmVec.SerializeToString(&dataString);

            int pos = dataString.size();

            // copy to buffer
            char tmpbuf[2048];
            std::strcpy(tmpbuf, dataString.c_str());
            memcpy (buf, tmpbuf, pos);
            buf_len = pos;
            
            return pos;
        }

        // done
        int CrosslinkSpectralMatch::Pack(char *buf, size_t &buf_len, CrosslinkSpectralMatch *csm)
        {
            std::vector<CrosslinkSpectralMatch*> csmVec;
            csmVec.push_back(csm);
   
            auto betaPeptide = csm->getBetaPeptide();
            if ( betaPeptide != nullptr ) {
                csvVec.push_back(betaPeptide);               
            }
            
            int pos = CrosslinkSpectralMatch::Pack(buf, buf_len, csmVec);

            return pos;            
        }
        
        // done
        SerializedCrosslinkSpectralMatch CrosslinkSpectralMatch::Pack_internal(CrosslinkSpectralMatch *csm)
        {
            // serialized csm object we will build
            SerializedCrosslinkSpectralMatch sCsm;

            auto mFrIons = csm->getMatchedFragmentIons ();
            auto dp = csm->digestionParams;
            auto uMapPep = csm->getPeptidesToMatchingFragments();
            std::vector<int> lPositions  = csm->getLinkPositions();
            std::vector<int> xlRanks = csm->getXlRank();
            bool has_beta_peptide = csm->getBetaPeptide() != nullptr;       

            sCsm.set_hasnotchvalue(csm->getNotch().has_value());   
            if ( csm->getNotch().has_value() ) {
                sCsm.set_notchvalue(csm->getNotch().value());
            }
            
            sCsm.set_xltotaltcore(csm->getXLTotalScore());
            sCsm.set_deltascore(sm->getDeltaScore());
            sCsm.set_score(csm->getScore());
            sCsm.set_runnerupscore(csm->getRunnerUpScore());
            sCsm.set_peptidemonoisotopicmass(csm->getPeptideMonisotopicMass().value());
            
            sCsm.set_scannumber(csm->getScanNumber());
            sCsm.set_xlproteinpos(csm->getXlProteinPos());
            sCsm.set_matchedfragmentionssize((int)csm->getMatchedFragmentIons().size());
            sCsm.set_lpositionssize((int)lPositions.size());
            sCsm.set_xlranksrize((int)xlRanks.size());

            sCsm.set_hasbetapeptide(has_beta_peptide);
                        
            PsmCrossType ctype = csm->getCrossType();
            sCsm.set_psmcrosstypeasstring(PsmCrossTypeToString(ctype));

            //Information required to replace the Scan datastructure
            auto tvar = csm->getPrecursorScanNumber();
            sCsm.set_hasprecursorscannumber(tvar.has_value());
            if (tvar.has_value()) {
                sCsm.set_precursorscannumber(tvar.value());
            }

            sCsm.set_scanexperimentalpeaks(csm->getScanExperimentalPeaks());
            sCsm.set_scanprecursorcharge(csm->getScanPrecursorCharge());

            sCsm.set_scanretentiontime(csm->getScanRetentionTime());
            sCsm.set_totalioncurrent(csm->getTotalIonCurrent());
            sCsm.set_scanprecursormonoisotopicpeakmz(csm->getScanPrecursorMonoisotopicPeakMz());
            sCsm.set_scanprecursormass(csm->getScanPrecursorMass());

            sCsm.set_fullfilepath(csm->getFullFilePath());

            FdrInfo *fdr = csm->getFdrInfo();

            // this routine sets all the required aspects of a packed line (e.g. header, length)
            SerializedFdrInfo sFdr = FdrInfo::Pack(fdr);
            sCsm.set_fdrinfo(sFdr);

            for (auto lp: lPositions) {
                sCsm.add_lpositions() = lp;
            }

            for ( auto xl: xlRanks) {
                sCsm.add_xlranks() = xl;
            }
            
            std::string s = dp->ToString();
            sCsm.set_digestionparamsstring(s);
                        
            //line 6-10: PeptideWithSetModifications;
            //Assuming right now only a single PeptideWithSetModifications
            if ( uMapPep.size() != 1 ) {
                std::cout << "CrosslinkSpectralMatch::Pack: Error - unordered_map has more than one entry!\n";
            }
            auto pep = std::get<0>(*uMapPep.begin());

            SerializedPeptideWithSetModifications sPep = PeptideWithSetModifications::Pack(pep);
            sCsm.set_peptide(sPep);
            
            //line 11-x: one line for each MatchedFragmentIon
            for ( auto i=0; i< mFrIons.size(); i++ ) {
                SerializedMatchedFragmentIon sMaF = MatchedFragmentIon::Pack(mFrIons[i]);
                sCsm.add_ions() = sMaF;
            }
            
            return sCsm;
        }

        /**
         * TODO...
         * 
         * Need to create an additional Pack() and Unpack()
         * method for each of the 3 objects packed/unpacked 
         * in this class. They each must be of the form:
         * 
         *      SerializedObjectType Pack(objectType obj);
         * 
         *      and
         * 
         *      void Unpack(objectType obj, pointerType destinationPointer);
        */

        // done
        void CrosslinkSpectralMatch::Unpack (char *buf, size_t buf_len, int count, size_t &len,
                                             std::vector<CrosslinkSpectralMatch*> &pepVec,
                                             const std::vector<Modification*> &mods,
                                             const std::vector<Protein *> &proteinList )
        {
            // convert input buf to string
            std::string dataString(buf);
            
            // parse string to object vector
            SerializedCrosslinkSpectralMatchVec sCsmVec;
            sCsmVec.ParseFromString(dataString);

            // convert object vec to CSM vector
            std::vector<SerializedCrosslinkSpectralMatch> tempCsmVec = {sCsmVec.csms().begin(), sCsmVec.csms().end()};

            int csmCount = tempCsmVec.size();
            int index = 0;

            // loop through CSM vector
            while (index < csmCount) {

                // unpack CSM
                CrosslinkSpectralMatch *pep;
                bool has_beta_peptide=false;
                CrosslinkSpectralMatch::Unpack_internal(tempCsmVec[index], &pep, mods, proteinList, has_beta_peptide);
                pepVec.push_back(pep);
                index++;

                // if CSM has betapep
                if (has_beta_peptide) {

                    // unpack betapep, set as CSM's beta peptide
                    CrosslinkSpectralMatch *beta_pep;
                    CrosslinkSpectralMatch::Unpack_internal(tempCsmVec[index], &pep, mods, proteinList, has_beta_peptide);
                    pep->setBetaPeptide(beta_pep);
                    index++;
                }
            }
        }

        // done
        void CrosslinkSpectralMatch::Unpack (char *buf, size_t buf_len, size_t &len,
                                             CrosslinkSpectralMatch** newCsm,
                                             const std::vector<Modification*> &mods,
                                             const std::vector<Protein *> &proteinList )
        {           
            // create CSM vector, call Unpack's vector version
            std::vector<CrosslinkSpectralMatch*> csmVec;
            CrosslinkSpectralMatch::Unpack(buf, buf_len, len, csmVec, mods, proteinList);
            if (csmVec.size() > 1) {
                *newCsm = csmVec[0];
            }
        }

        // done
        void CrosslinkSpectralMatch::Unpack_internal (SerializedCrosslinkSpectralMatch sCsm,
                                                      CrosslinkSpectralMatch** newCsm,
                                                      const std::vector<Modification*> &mods,
                                                      const std::vector<Protein *> &proteinList,
                                                      bool &has_beta_peptide )
        {            
            int notch=-1, scannumber, proteinPos, matchedFragmentIonsVecsize, lpositionsize, xlranksize;
            double  deltaScore, XLTotalScore, score, runnerUpScore, peptideMonisotopicMass;
            bool  has_notch;

            has_notch = sCsm.hasnotchvalue();
            if (has_notch) {
                notch = sCsm.notchvalue();
            }

            XLTotalScore = sCsm.xltotalscore();
            deltaScore = sCsm.seltascore();
            score = sCsm.score();
            runnerUpScore =sCsm.runnerupscore();
            peptideMonisotopicMass = sCsm.peptidemonoisotopicmass();
            
            scannumber = sCsm.scannumber();
            proteinPos = sCsm.xlproteinpos();
            matchedFragmentIonsVecsize = sCsm.matchedfragmentionssize();
            lpositionsize = sCsm.lpositionssize();
            xlranksize = sCsm.xlrankssize();
            
            has_beta_peptide = sCsm.hasbetapeptide();

            std::string tmpstring = sCsm.psmcrosstypeasstring();
            PsmCrossType ctype = PsmCrossTypeFromString(tmpstring);

            //Information required to replace the Scan datastructure
            bool has_tvar;
            int  itvar;
            std::optional<int> scanPrecursorScanNumber;

            has_tvar = sCsm.hasprecursorscannumber();
            if (has_tvar) {
                itvar = sCsm.precursorscannumber();
                scanPrecursorScanNumber = std::make_optional(itvar);
            }

            int scanExperimentalPeaks = sCsm.scanexperimentalpeaks();
            int scanPrecursorCharge = sCsm.scanprecursorcharge();

            double scanRetentionTime = sCsm.scanretentiontime();
            double scanTotalIonCurrent = sCsm.totalioncurrent();
            double scanPrecursorMonoisotopicPeakMz = sCsm.scanprecursormonoisotopicpeakmz();
            double scanPrecursorMass = sCsm.scanprecursormass();

            std::string scanFullFilePath = sCsm.fullfilepath();;

            FdrInfo* fdr = nullptr;
            FdrInfo::Unpack(sCsm.fdrinfo(), &fdr);

            std::vector<int> linkPosvec = {sCsm.lpositions().begin(), sCsm.lpositions().end()};           
            std::vector<int> xlRankVec = {sCsm.xranks().begin(), sCsm.xranks().end()};

            std::string dpstring = sCsm.digestionparamsstring();
            DigestionParams *dp = DigestionParams::FromString(dpstring);           

            //line 6-10: PeptideWithSetModifications
            PeptideWithSetModifications* pep;
            PeptideWithSetModifications::Unpack(sCsm.peptide(), &pep);
            pep->SetNonSerializedPeptideInfo(mods, proteinList);

#ifdef DEBUG
            // Safety check:
            double monIsotopicMass = pep->getMonoisotopicMass();
            if (monIsotopicMass != peptideMonisotopicMass) {
                std::cout << "Safety check failed when reconstructing PeptideWithSetMOdifications in CrosslinkSpectralMatch::Unpack(). " <<
                    "monIsotopicMass is " << monIsotopicMass << " should be " << peptideMonisotopicMass << std::endl;
            }
#endif
            
            // get serialized MaFIons from serialized CSM
            std::vector<SerializedMatchedFragmention> tempMaFVec = {sCsm.ions().begin(), sCsm.ions().end()};

            // declare new MaFIon vec
            std::vector<MatchedFragmentIon*> matchedFragmentIonsVec;

            // unpack each serialized MaF and push to MaF vector
            for (auto i = 0; i < matchedFragmentIonsVecsize; i++) {
                MatchedFragmentIon *ion;
                MatchedFragmentIon::Unpack(tempMaFVec[i], &ion);
                matchedFragmentIonsVec.push_back(ion);                    
            }

            // We are treating scannumber and scanindex as the same here. First, it is actually really the same
            // in many scenarios. Second, scanindex is not really used for the subsequent operations as far
            // as I can see right now.
            CrosslinkSpectralMatch *csm = new CrosslinkSpectralMatch ( pep, notch, XLTotalScore, scannumber,
                                                                       scanFullFilePath, scannumber,
                                                                       scanPrecursorScanNumber,
                                                                       scanRetentionTime, scanExperimentalPeaks,
                                                                       scanTotalIonCurrent, scanPrecursorCharge,
                                                                       scanPrecursorMonoisotopicPeakMz,
                                                                       scanPrecursorMass, dp,
                                                                       matchedFragmentIonsVec );
            csm->setXLTotalScore(XLTotalScore);
            csm->setDeltaScore(deltaScore);
            csm->setXlProteinPos(proteinPos);
            csm->setScore(score);
            csm->setRunnerUpScore(runnerUpScore);
            csm->setCrossType (ctype);
            csm->setXlRank(xlRankVec);
            csm->setLinkPositions(linkPosvec);
            if (fdr != nullptr) {
                csm->setFdrInfo(fdr);
            }
            csm->ResolveAllAmbiguities();

            // This is a hack. Need to find a proper solution.
            // In some situations involving Mods, the PeptideMonoisotopicMass
            // is not correct after deserialization.
            csm->setPeptideMonisotopicMass(std::make_optional(peptideMonisotopicMass));
            
            *newCsm = csm;
            
            return;
        }   
    }
}
