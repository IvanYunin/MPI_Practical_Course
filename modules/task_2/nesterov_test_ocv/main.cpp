// Copyright [year] <Copyright Owner>
#include <mpi.h>
#include <assert.h>
#include <opencv2/opencv.hpp>
#include <iostream>

int main(int argc, char **argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    cv::Mat  input(cv::Size(1920, 1080), CV_8UC1, cv::Scalar(0));
    cv::Mat output(cv::Size(1920, 1080), CV_8UC1, cv::Scalar(0));
    cv::RNG rng;
    rng.fill(input, cv::RNG::UNIFORM, cv::Scalar(0), cv::Scalar(255));
    cv::adaptiveThreshold(input, output, 255, cv::ADAPTIVE_THRESH_MEAN_C,
                          cv::THRESH_BINARY, 3, 1.5);
    std::cout << "Count non zero element = ";
    std::cout << cv::countNonZero(output) << '\n';

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
