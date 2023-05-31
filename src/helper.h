#pragma once
#ifndef __HELPER_H__
#define __HELPER_H__

#include"distmesh.h"
#include <fstream>
#include <chrono>

namespace distmesh {
    namespace helper {
        // save eigen array to text file
        template <typename type>
        void savetxt(Eigen::Ref<Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> const> const array,
            std::string const& filename) {
            // open file
            std::ofstream file(filename);

            // save array to file with high precision
            Eigen::IOFormat const format(Eigen::FullPrecision, Eigen::DontAlignCols);
            file << array.format(format);

            file.close();
        }

        class HighPrecisionTime {
        private:
            std::chrono::high_resolution_clock::time_point time;

        public:
            HighPrecisionTime() {
                this->restart();
            }

            void restart() {
                this->time = std::chrono::high_resolution_clock::now();
            }

            double elapsed() const {
                return std::chrono::duration_cast<std::chrono::duration<double>>(
                    std::chrono::high_resolution_clock::now() - this->time).count();
            }
        };
    }
}

#endif