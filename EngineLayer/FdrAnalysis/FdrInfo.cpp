#include "FdrInfo.h"
#include "FdrInfo_generated.h"

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
	    flatbuffers::Offset<SerializedFdrInfo> FdrInfo::Pack(FdrInfo *fdr)
	    {
		flatbuffers::FlatBufferBuilder builder;
                bool hasFdr = (fdr != nullptr);

                // build serialized fdr if it exists
                if (hasFdr) {
                    double cumulativeTarget = fdr->getCumulativeTarget();
                    double cumulativeDecoy = fdr->getCumulativeDecoy();
                    double qValue = fdr->getQValue();
                    double cumulativeTargetNotch = fdr->getCumulativeTargetNotch();
                    double cumulativeDecoyNotch = fdr->getCumulativeDecoyNotch();
                    double qValueNotch = fdr->getQValueNotch();;
                    double maximumLikelihood = fdr->getMaximumLikelihood();
                    double eValue = fdr->getEValue();
                    double eScore = fdr->getEScore();
                    bool calculateEValue = fdr->getCalculateEValue();

                    auto fdrInfo = CreateSerializedFdrInfo(builder, cumulativeTarget, cumulativeDecoy, qValue,
                                                        cumulativeTargetNotch, cumulativeDecoyNotch,
                                                        qValueNotch, maximumLikelihood, eValue, eScore,
                                                        calculateEValue, hasFdr);
                    return fdrInfo;
                }
                else {
                    auto fdrInfo = CreateSerializedFdrInfo(builder, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false, hasFdr);
                    return fdrInfo;
                }
	    }
            
            int FdrInfo::Pack( char* buf, size_t &buf_len, FdrInfo *fdr)
            {
                flatbuffers::FlatBufferBuilder builder;
                bool hasFdr = (fdr != nullptr);
                
                // build serialized fdr if it exists
                if (hasFdr) {
                    double cumulativeTarget = fdr->getCumulativeTarget();
                    double cumulativeDecoy = fdr->getCumulativeDecoy();
                    double qValue = fdr->getQValue();
                    double cumulativeTargetNotch = fdr->getCumulativeTargetNotch();
                    double cumulativeDecoyNotch = fdr->getCumulativeDecoyNotch();
                    double qValueNotch = fdr->getQValueNotch();;
                    double maximumLikelihood = fdr->getMaximumLikelihood();
                    double eValue = fdr->getEValue();
                    double eScore = fdr->getEScore();
                    bool calculateEValue = fdr->getCalculateEValue();

                    auto fdrInfo = CreateSerializedFdrInfo(builder, cumulativeTarget, cumulativeDecoy, qValue,
                                                        cumulativeTargetNotch, cumulativeDecoyNotch,
                                                        qValueNotch, maximumLikelihood, eValue, eScore,
                                                        calculateEValue, hasFdr);
		    builder.Finish(fdrInfo);
                }
                else {
                    auto fdrInfo = CreateSerializedFdrInfo(builder, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false, hasFdr);
                    builder.Finish(fdrInfo);
		}
                
                char* tempBuf = (char*)builder.GetBufferPointer();

                int pos = builder.GetSize();
                
                // copy to buffer
                memcpy(buf, tempBuf, pos);
                buf_len = pos;
                return pos;
            }

            // for use in CSM unpacking
            void FdrInfo::Unpack(const SerializedFdrInfo* sFdrInfo, FdrInfo **newfdr)
            {
                double cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch;
                double cumulativeDecoyNotch, qValueNotch, maximumLikelihood, eValue, eScore;
                bool calculateEValue = false;
                bool has_fdr = false;

               	has_fdr = sFdrInfo->hasFdr();

                // rebuild fdr
                if (has_fdr) {
                    FdrInfo* tempVar = new FdrInfo();
                    tempVar->setCumulativeTarget(sFdrInfo->cumulativeTarget());
                    tempVar->setCumulativeDecoy(sFdrInfo->cumulativeDecoy());
                    tempVar->setQValue(sFdrInfo->qValue());
                    tempVar->setCumulativeTargetNotch(sFdrInfo->cumulativeTargetNotch());
                    tempVar->setCumulativeDecoyNotch(sFdrInfo->cumulativeDecoyNotch());
                    tempVar->setQValueNotch(sFdrInfo->qValueNotch());
                    tempVar->setMaximumLikelihood(sFdrInfo->maximumLikelihood());
                    tempVar->setEScore(sFdrInfo->eScore());
                    tempVar->setEValue(sFdrInfo->eValue());
                    tempVar->setCalculateEValue(sFdrInfo->calculateEValue());

                    *newfdr = tempVar;
                }
                else {
                    *newfdr = nullptr;
                }
            }

            void FdrInfo::Unpack(char* buf, size_t &len, FdrInfo **newfdr)
            {
                auto sFdrInfo = GetSerializedFdrInfo((uint8_t*)buf);        
                bool has_fdr = sFdrInfo->hasFdr();

                // rebuild fdr
                if (has_fdr) {
                    FdrInfo* tempVar= new FdrInfo();
                    tempVar->setCumulativeTarget(sFdrInfo->cumulativeTarget());
                    tempVar->setCumulativeDecoy(sFdrInfo->cumulativeDecoy());
                    tempVar->setQValue(sFdrInfo->qValue());
                    tempVar->setCumulativeTargetNotch(sFdrInfo->cumulativeTargetNotch());
                    tempVar->setCumulativeDecoyNotch(sFdrInfo->cumulativeDecoyNotch());
                    tempVar->setQValueNotch(sFdrInfo->qValueNotch());
                    tempVar->setMaximumLikelihood(sFdrInfo->maximumLikelihood());
                    tempVar->setEScore(sFdrInfo->eScore());
                    tempVar->setEValue(sFdrInfo->eValue());
                    tempVar->setCalculateEValue(sFdrInfo->calculateEValue());

                    *newfdr = tempVar;
                }
                else {
                    *newfdr = nullptr;
                }
            }
        }
}
