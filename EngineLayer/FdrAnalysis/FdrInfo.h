#pragma once

#include <string>
#include <msgpack.hpp>

namespace EngineLayer
{
    namespace FdrAnalysis
    {
	struct SerializedFdrInfo {

            double cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch;
            double cumulativeDecoyNotch, qValueNotch, maximumLikelihood, eValue, eScore;
            bool calculateEValue, has_fdr;

            MSGPACK_DEFINE(cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch,
                            cumulativeDecoyNotch, qValueNotch, maximumLikelihood, eValue, eScore,
                            calculateEValue, has_fdr);
        };

        class FdrInfo
        {
        private:
            double privateCumulativeTarget = 0;
            double privateCumulativeDecoy = 0;
            double privateCumulativeTargetNotch = 0;
            double privateCumulativeDecoyNotch = 0;
            double privateQValue = 0;
            double privateQValueNotch = 0;
            bool privateCalculateEValue = false;
            double privateMaximumLikelihood = 0;
            double privateEValue = 0;
            double privateEScore = 0;
            
        public:
            double getCumulativeTarget() const;
            void setCumulativeTarget(double value);
            double getCumulativeDecoy() const;
            void setCumulativeDecoy(double value);
            double getCumulativeTargetNotch() const;
            void setCumulativeTargetNotch(double value);
            double getCumulativeDecoyNotch() const;
            void setCumulativeDecoyNotch(double value);
            double getQValue() const;
            void setQValue(double value);
            double getQValueNotch() const;
            void setQValueNotch(double value);
            bool getCalculateEValue() const;
            void setCalculateEValue(bool value);
            double getMaximumLikelihood() const;
            void setMaximumLikelihood(double value);
            double getEValue() const;
            void setEValue(double value);
            double getEScore() const;
            void setEScore(double value);

            /// <summary>
            /// Pack an FdrInfo into a character buffer.
            /// Required for Communication Operations in HPCMetaMorpheus
            ///
            /// Arguments:
            /// buf :     INOUT buffer used for packing
            /// buf_size: IN size of the allocated buffer provided by the upper layer
            ///           OUT size of required buffer if not large enough (return value -1)
            ///               or number of bytes used for packgin (return value > 0)
            /// fdr :     IN pointer to FdrInfo to pack. Can be a nullptr
            ///
            /// Return value:
            ///   -1 : input buffer was not large enough. buf_size will contain the required number
            ///        of bytes in this case
            ///   >0 : packing successful, number of bytes used up.
            /// </summary>
            static const char* Pack(FdrInfo *fdr); 

            /// <summary>
            /// Functionality used to reconstruct an FdrInfo based on a
            /// packed buffer.
            ///
            /// Arguments
            /// ---------
            /// buf:      IN input character buffer
            /// len:      OUT number of bytes used for unpacking an FdrInfo
            /// newFdr:   OUT new FdrInfo(s). Can be nullptr.
            ///
            /// </summary>
            static void Unpack(char* buf, FdrInfo **newfdr);
        };
    }
}
