﻿#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <tuple>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }

using namespace Proteomics;
using namespace Proteomics::AminoAcidPolymer;

namespace EngineLayer
{
	namespace HistogramAnalysis
	{
		class Bin
		{
		private:
			int privatePepNlocCount = 0;
			int privatePepClocCount = 0;
			int privateProtNlocCount = 0;
			int privateProtClocCount = 0;
			std::string privateCombos = "-";
			std::string privateUnimodDiffs = "-";
			std::string privateUniprotID = "-";
			std::string privateUnimodFormulas = "-";
			std::string privateUnimodId = "-";
			double privateMassShift = 0;
			std::string privateMine;
			std::unordered_map<char, int> privateAAsInCommon;
			int privateOverlapping = 0;
			double privateFracWithSingle = 0;
			double privateMedianLength = 0;

		public:
			std::string AA = "-";
			std::unordered_map<char, int> ResidueCount;
			std::unordered_map<std::string, std::tuple<std::string, std::string, PeptideSpectralMatch*>> UniquePSMs;
			std::unordered_map<std::string, int> ModsInCommon;

			Bin(double massShift);

				int getPepNlocCount() const;
				void setPepNlocCount(int value);
				int getPepClocCount() const;
				void setPepClocCount(int value);
				int getProtNlocCount() const;
				void setProtNlocCount(int value);
				int getProtClocCount() const;
				void setProtClocCount(int value);
				std::string getCombos() const;
				void setCombos(const std::string &value);
				std::string getUnimodDiffs() const;
				void setUnimodDiffs(const std::string &value);
				std::string getUniprotID() const;
				void setUniprotID(const std::string &value);
				std::string getUnimodFormulas() const;
				void setUnimodFormulas(const std::string &value);
				std::string getUnimodId() const;
				void setUnimodId(const std::string &value);
				double getMassShift() const;

			int getCount() const;

			int getCountDecoy() const;

			int getCountTarget() const;

			int getLocalizeableTarget() const;

				std::string getMine() const;
				void setMine(const std::string &value);
				std::unordered_map<char, int> getAAsInCommon() const;
				void setAAsInCommon(const std::unordered_map<char, int> &value);
				int getOverlapping() const;
				void setOverlapping(int value);
				double getFracWithSingle() const;
				void setFracWithSingle(double value);
				double getMedianLength() const;
				void setMedianLength(double value);

			void IdentifyResidues();

			void IdentifyCombos(double v, std::unordered_set<std::tuple<double, double, double>> &ok);

			double ComputeZ(double v);

			void IdentifyUniprotBins(double v);

			void IdentifyAA(double v);

			void IdentifyUnimodBins(double v);

			void Add(PeptideSpectralMatch *ok);
		};
	}
}
