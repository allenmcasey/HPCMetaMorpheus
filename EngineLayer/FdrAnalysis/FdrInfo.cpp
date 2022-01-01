#include "FdrInfo.h"

#include <string.h>
#include <sstream>
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
            
            const char* FdrInfo::Pack(FdrInfo *fdr)
            {
		std::stringstream sstream;
		bool has_fdr = (fdr != nullptr);
    								
    		if (has_fdr) {
    	
    			SerializedFdrInfo sFdr = {
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
    			msgpack::pack(sstream, sFdr);
    			const std::string& tmp = sstream.str();
			const char* buf = tmp.c_str();
			return buf;
    		}
    		return nullptr;
            }
	
	    void FdrInfo::Unpack(char* buf, FdrInfo **newfdr)
      	    {
            	double cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch;
            	double cumulativeDecoyNotch, qValueNotch, maximumLikelihood, eValue, eScore;
            	bool calculateEValue=false;
            	bool has_fdr = false;

            	std::stringstream buffer;
            	buffer << buf;

            	msgpack::object_handle objHandle = msgpack::unpack(buffer.str().data(), buffer.str().size());
            	msgpack::object const& obj = objHandle.get();
            	auto sFdrInfo = obj.as<SerializedFdrInfo>();

            	has_fdr = sFdrInfo.has_fdr;

            	if (has_fdr) {

                	FdrInfo* tempVar= new FdrInfo();

                	tempVar->setCumulativeTarget(sFdrInfo.cumulativeTarget);
                	tempVar->setCumulativeDecoy(sFdrInfo.cumulativeDecoy);
                	tempVar->setQValue(sFdrInfo.qValue);
                	tempVar->setCumulativeTargetNotch(sFdrInfo.cumulativeTargetNotch);
                	tempVar->setCumulativeDecoyNotch(sFdrInfo.cumulativeDecoyNotch);
                	tempVar->setQValueNotch(sFdrInfo.qValueNotch);
                	tempVar->setMaximumLikelihood(sFdrInfo.maximumLikelihood);
                	tempVar->setEScore(sFdrInfo.eScore);
                	tempVar->setEValue(sFdrInfo.eValue);
                	tempVar->setCalculateEValue(sFdrInfo.calculateEValue);

                	*newfdr = tempVar;
            	}
            	else {
                	*newfdr = nullptr;
            	}
            }
        }
}
