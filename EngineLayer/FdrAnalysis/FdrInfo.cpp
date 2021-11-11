#include "FdrInfo.pb.h"
#include "FdrInfo.h"

#include <string.h>
#include "BinaryPack.h"

namespace EngineLayer
{
	namespace FdrAnalysis
	{

            double FdrInfo::getCumulativeTarget() const
            {
                return privateCumulativeTarget;
            }
            
            void FdrInfo::setCumulativeTarget(double value)
            {
                privateCumulativeTarget = value;
            }
            
            double FdrInfo::getCumulativeDecoy() const
            {
                return privateCumulativeDecoy;
            }
            
            void FdrInfo::setCumulativeDecoy(double value)
            {
                privateCumulativeDecoy = value;
            }
            
            double FdrInfo::getCumulativeTargetNotch() const
            {
                return privateCumulativeTargetNotch;
            }
            
            void FdrInfo::setCumulativeTargetNotch(double value)
            {
                privateCumulativeTargetNotch = value;
            }
            
            double FdrInfo::getCumulativeDecoyNotch() const
            {
                return privateCumulativeDecoyNotch;
            }
            
            void FdrInfo::setCumulativeDecoyNotch(double value)
            {
                privateCumulativeDecoyNotch = value;
            }
            
            double FdrInfo::getQValue() const
            {
                return privateQValue;
            }
            
            void FdrInfo::setQValue(double value)
            {
                privateQValue = value;
            }
            
            double FdrInfo::getQValueNotch() const
            {
                return privateQValueNotch;
            }
            
            void FdrInfo::setQValueNotch(double value)
            {
                privateQValueNotch = value;
            }
            
            bool FdrInfo::getCalculateEValue() const
            {
                return privateCalculateEValue;
            }
            
            void FdrInfo::setCalculateEValue(bool value)
            {
                privateCalculateEValue = value;
            }
            
            double FdrInfo::getMaximumLikelihood() const
            {
                return privateMaximumLikelihood;
            }
            
            void FdrInfo::setMaximumLikelihood(double value)
            {
                privateMaximumLikelihood = value;
            }
            
            double FdrInfo::getEValue() const
            {
                return privateEValue;
            }
            
            void FdrInfo::setEValue(double value)
            {
                privateEValue = value;
            }
            
            double FdrInfo::getEScore() const
            {
                return privateEScore;
            }
            
            void FdrInfo::setEScore(double value)
            {
                privateEScore = value;
            }

            // for use in CSM packing
            SerializedFdrInfo FdrInfo::Pack(FdrInfo *fdr)
            {
                SerializedFdrInfo sFdrInfo;

                sFdrInfo.set_hasFdr(fdr != nullptr);
                
                // build serialized fdr if it exists
                if (fdr != nullptr) {
                    sFdrInfo.set_cumulativeTarget(fdr->getCumulativeTarget());
                    sFdrInfo.set_cumulativeDecoy(fdr->getCumulativeDecoy());
                    sFdrInfo.set_qValue(fdr->getQValue());
                    sFdrInfo.set_cumulativeTargetNotch(fdr->getCumulativeTargetNotch());
                    sFdrInfo.set_cumulativeDecoyNotch(fdr->getCumulativeDecoyNotch());
                    sFdrInfo.set_qValueNotch(fdr->getQValueNotch());
                    sFdrInfo.set_maximumLikelihood(fdr->getMaximumLikelihood());
                    sFdrInfo.set_eValue(fdr->getEValue());
                    sFdrInfo.set_eScore(fdr->getEScore());
                    sFdrInfo.set_calculateEValue(fdr->getCalculateEValue());
                }

                return sFdrInfo;
            }
            
            int FdrInfo::Pack(char* buf, size_t &buf_len, FdrInfo *fdr)
            {
                SerializedFdrInfo sFdrInfo;

                sFdrInfo.set_hasFdr(fdr != nullptr);
                
                // build serialized fdr if it exists
                if (fdr != nullptr) {
                    sFdrInfo.set_cumulativeTarget(fdr->getCumulativeTarget());
                    sFdrInfo.set_cumulativeDecoy(fdr->getCumulativeDecoy());
                    sFdrInfo.set_qValue(fdr->getQValue());
                    sFdrInfo.set_cumulativeTargetNotch(fdr->getCumulativeTargetNotch());
                    sFdrInfo.set_cumulativeDecoyNotch(fdr->getCumulativeDecoyNotch());
                    sFdrInfo.set_qValueNotch(fdr->getQValueNotch());
                    sFdrInfo.set_maximumLikelihood(fdr->getMaximumLikelihood());
                    sFdrInfo.set_eValue(fdr->getEValue());
                    sFdrInfo.set_eScore(fdr->getEScore());
                    sFdrInfo.set_calculateEValue(fdr->getCalculateEValue());
                }
                
                // serialize FdrInfo
                std::string dataString;
                sFdrInfo.SerializeToString(&dataString);
                int pos = dataString.size();
                
                // copy to buffer
                char tmpbuf[128];
                memset(tmpbuf, 0, 128);
                std::strcpy(tmpbuf, dataString.c_str());
                memcpy(buf, tmpbuf, pos);
                buf_len = pos;
                return pos;
            }

            // for use in CSM unpacking
            void FdrInfo::Unpack(SerializedFdrInfo sFdrInfo, FdrInfo **newfdr)
            {
                double cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch;
                double cumulativeDecoyNotch, qValueNotch, maximumLikelihood, eValue, eScore;
                bool calculateEValue = false;
                bool has_fdr = false;

                // check if the fdr is null
                has_fdr = sFdrInfo.hasfdr();
                
                // rebuild fdr
                if (has_fdr) {
                    FdrInfo* tempVar= new FdrInfo();
                    tempVar->setCumulativeTarget(sFdrInfo.cumulativetarget());
                    tempVar->setCumulativeDecoy(sFdrInfo.cumulativedecoy());
                    tempVar->setQValue(sFdrInfo.qvalue());
                    tempVar->setCumulativeTargetNotch(sFdrInfo.cumulativetargetnotch());
                    tempVar->setCumulativeDecoyNotch(sFdrInfo.cumulativedecoynotch());
                    tempVar->setQValueNotch(sFdrInfo.qvaluenotch());
                    tempVar->setMaximumLikelihood(sFdrInfo.maximumlikelihood());
                    tempVar->setEScore(sFdrInfo.escore());
                    tempVar->setEValue(sFdrInfo.evalue());
                    tempVar->setCalculateEValue(sFdrInfo.calculateevalue());

                    *newfdr = tempVar;
                }
                else {
                    *newfdr = nullptr;
                }
            }

            void FdrInfo::Unpack(char* buf, size_t &len, FdrInfo **newfdr)
            {
                double cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch;
                double cumulativeDecoyNotch, qValueNotch, maximumLikelihood, eValue, eScore;
                bool calculateEValue = false;
                bool has_fdr = false;
                
                // convert input buf to string
                std::string dataString(buf);
                len = dataString.size();

                // parse fdr object from string
                SerializedFdrInfo sFdrInfo;
                sFdrInfo.ParseFromString(dataString);

                // check if the fdr is null
                has_fdr = sFdrInfo.hasfdr();
                
                // rebuild fdr
                if (has_fdr) {
                    FdrInfo* tempVar= new FdrInfo();
                    tempVar->setCumulativeTarget(sFdrInfo.cumulativetarget());
                    tempVar->setCumulativeDecoy(sFdrInfo.cumulativedecoy());
                    tempVar->setQValue(sFdrInfo.qvalue());
                    tempVar->setCumulativeTargetNotch(sFdrInfo.cumulativetargetnotch());
                    tempVar->setCumulativeDecoyNotch(sFdrInfo.cumulativedecoynotch());
                    tempVar->setQValueNotch(sFdrInfo.qvaluenotch());
                    tempVar->setMaximumLikelihood(sFdrInfo.maximumlikelihood());
                    tempVar->setEScore(sFdrInfo.escore());
                    tempVar->setEValue(sFdrInfo.evalue());
                    tempVar->setCalculateEValue(sFdrInfo.calculateevalue());

                    *newfdr = tempVar;
                }
                else {
                    *newfdr = nullptr;
                }
            }
        }
}