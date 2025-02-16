
#include <sweet/Memory/parmemcpy.hpp>

#include "DataGrid.hpp"

namespace sweet {
    namespace Data {
        namespace Sphere2D {


            void DataGrid::swap(
                    DataGrid &i_sphere2DData
            ) {
                SWEET_ASSERT(sphere2DDataConfig == i_sphere2DData.sphere2DDataConfig);

                std::swap(grid_space_data, i_sphere2DData.grid_space_data);
            }


            DataGrid::DataGrid(
                    const Sphere2D::Config &i_sphere2DDataConfig
            ) :
                    DataGrid(&i_sphere2DDataConfig) {
            }


            DataGrid::DataGrid(
                    const Sphere2D::Config *i_sphere2DDataConfig
            ) :
            //! important: set this to nullptr, since a check for this will be performed by setup(...)
                    sphere2DDataConfig(i_sphere2DDataConfig),
                    grid_space_data(nullptr) {
                alloc_data();
            }


            DataGrid::DataGrid(
                    const Sphere2D::Config *i_sphere2DDataConfig,
                    double i_value
            ) :
            //! important: set this to nullptr, since a check for this will be performed by setup(...)
                    sphere2DDataConfig(i_sphere2DDataConfig),
                    grid_space_data(nullptr) {
                alloc_data();
                grid_setValue(i_value);
            }


            DataGrid::DataGrid() :
                    sphere2DDataConfig(nullptr),
                    grid_space_data(nullptr) {
            }


            DataGrid::DataGrid(
                    const DataGrid &i_sph_data
            ) :
                    sphere2DDataConfig(i_sph_data.sphere2DDataConfig),
                    grid_space_data(nullptr) {
                if (i_sph_data.sphere2DDataConfig != nullptr)
                    alloc_data();

                operator=(i_sph_data);
            }


            DataGrid::DataGrid(
                    DataGrid &&i_sph_data
            ) :
                    sphere2DDataConfig(i_sph_data.sphere2DDataConfig),
                    grid_space_data(nullptr) {
                if (i_sph_data.sphere2DDataConfig == nullptr)
                    return;

                std::swap(grid_space_data, i_sph_data.grid_space_data);
            }


/*!
 * Run validation checks to make sure that the physical and spectral spaces match in size
 */
            DataGrid::~DataGrid() {
                clear();
            }


            void DataGrid::clear() {
                if (grid_space_data != nullptr) {
                    sweet::Memory::MemBlockAlloc::free(grid_space_data,
                                                       sphere2DDataConfig->grid_number_elements * sizeof(double));
                    grid_space_data = nullptr;
                }

                sphere2DDataConfig = nullptr;
            }


            void DataGrid::check(
                    const Sphere2D::Config *i_sphere2DDataConfig
            ) const {
                SWEET_ASSERT(sphere2DDataConfig->grid_num_lat == i_sphere2DDataConfig->grid_num_lat);
                SWEET_ASSERT(sphere2DDataConfig->grid_num_lon == i_sphere2DDataConfig->grid_num_lon);
            }


            DataGrid &DataGrid::operator=(
                    const DataGrid &i_sph_data
            ) {
                if (i_sph_data.sphere2DDataConfig == nullptr)
                    return *this;

                if (sphere2DDataConfig == nullptr)
                    setup(i_sph_data.sphere2DDataConfig);

                sweet::Memory::parmemcpy(grid_space_data, i_sph_data.grid_space_data,
                                         sizeof(double) * sphere2DDataConfig->grid_number_elements);

                return *this;
            }


            DataGrid &DataGrid::operator=(
                    DataGrid &&i_sph_data
            ) {
                if (sphere2DDataConfig == nullptr)
                    setup(i_sph_data.sphere2DDataConfig);

                std::swap(grid_space_data, i_sph_data.grid_space_data);

                return *this;
            }


            DataGrid DataGrid::operator+(
                    double i_value
            ) const {
                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    out.grid_space_data[idx] = grid_space_data[idx] + i_value;

                return out;
            }


            DataGrid DataGrid::operator+(
                    const DataGrid &i_sph_data
            ) const {
                check(i_sph_data.sphere2DDataConfig);

                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    out.grid_space_data[idx] = grid_space_data[idx] + i_sph_data.grid_space_data[idx];

                return out;
            }


            DataGrid &DataGrid::operator+=(
                    const DataGrid &i_sph_data
            ) {
                check(i_sph_data.sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    grid_space_data[idx] += i_sph_data.grid_space_data[idx];

                return *this;
            }


            DataGrid &DataGrid::operator+=(
                    double i_scalar
            ) {
                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    grid_space_data[idx] += i_scalar;

                return *this;
            }


            DataGrid DataGrid::operator-(
                    double i_value
            ) const {
                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    out.grid_space_data[idx] = grid_space_data[idx] - i_value;

                return out;
            }


            DataGrid DataGrid::operator-() {
                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    out.grid_space_data[idx] = -grid_space_data[idx];

                return out;
            }


            DataGrid DataGrid::operator-(
                    const DataGrid &i_sph_data
            ) const {
                check(i_sph_data.sphere2DDataConfig);

                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    out.grid_space_data[idx] = grid_space_data[idx] - i_sph_data.grid_space_data[idx];

                return out;
            }


            DataGrid &DataGrid::operator-=(
                    double i_scalar
            ) {
                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    grid_space_data[idx] -= i_scalar;

                return *this;
            }


            DataGrid &DataGrid::operator-=(
                    const DataGrid &i_sph_data
            ) {
                check(i_sph_data.sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    grid_space_data[idx] -= i_sph_data.grid_space_data[idx];

                return *this;
            }


            DataGrid DataGrid::operator*(
                    const double i_value
            ) const {
                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
                    out.grid_space_data[i] = grid_space_data[i] * i_value;

                return out;
            }


            const DataGrid &DataGrid::operator*=(
                    const double i_value
            ) const {
                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    grid_space_data[idx] *= i_value;

                return *this;
            }


            DataGrid DataGrid::operator*(
                    const DataGrid &i_sph_data
            ) const {
                check(i_sph_data.sphere2DDataConfig);

                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
                    out.grid_space_data[i] = grid_space_data[i] * i_sph_data.grid_space_data[i];

                return out;
            }


            DataGrid DataGrid::operator/(
                    double i_value
            ) const {
                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    out.grid_space_data[idx] = grid_space_data[idx] / i_value;

                return out;
            }


            DataGrid DataGrid::operator/(
                    const DataGrid &i_sph_data
            ) const {
                check(i_sph_data.sphere2DDataConfig);

                check(i_sph_data.sphere2DDataConfig);

                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
                    out.grid_space_data[i] = grid_space_data[i] / i_sph_data.grid_space_data[i];

                return out;
            }


            DataGrid DataGrid::operator_scalar_sub_this(
                    double i_value
            ) const {
                DataGrid out(sphere2DDataConfig);

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    out.grid_space_data[idx] = i_value - grid_space_data[idx];

                return out;
            }


            bool DataGrid::setup(
                    const Sphere2D::Config *i_sphere2DDataConfig
            ) {
                if (sphere2DDataConfig != nullptr)
                    SWEETErrorFatal("Setup called twice!");

                sphere2DDataConfig = i_sphere2DDataConfig;
                return alloc_data();
            }


            bool DataGrid::setup(
                    const Sphere2D::Config &i_sphere2DDataConfig
            ) {
                return setup(&i_sphere2DDataConfig);
            }


            bool DataGrid::alloc_data() {
                SWEET_ASSERT(grid_space_data == nullptr);
                grid_space_data = sweet::Memory::MemBlockAlloc::alloc<double>(
                        sphere2DDataConfig->grid_number_elements * sizeof(double));
                return true;
            }


            void DataGrid::setup_if_required(
                    const Sphere2D::Config *i_sphere2DDataConfig
            ) {
                if (sphere2DDataConfig != nullptr)
                    if (sphere2DDataConfig == i_sphere2DDataConfig)
                        return;

                setup(i_sphere2DDataConfig);
            }


            void DataGrid::grid_update_lambda(
                    std::function<void(double, double, double &)> i_lambda
            ) {
                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD

#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS

                for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
                {
                    double lon_degree = ((double)i/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

                    for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
                    {
                        //double colatitude = acos(shtns->ct[j]);

                        /*
                         * Colatitude is 0 at the north pole and 180 at the south pole
                         *
                         * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
                         */
                        //double lat_degree = M_PI*0.5 - colatitude;
                        double lat_degree = sphere2DDataConfig->lat[j];

                        i_lambda(lon_degree, lat_degree, grid_space_data[i*sphere2DDataConfig->grid_num_lat + j]);
                    }
                }

#else

                for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++) {
                    double lat_degree = sphere2DDataConfig->lat[jlat];

                    for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++) {
                        double lon_degree = ((double) ilon / (double) sphere2DDataConfig->grid_num_lon) * 2.0 * M_PI;

                        //double colatitude = acos(shtns->ct[j]);

                        /*
                         * Colatitude is 0 at the north pole and 180 at the south pole
                         *
                         * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
                         */
                        //double lat_degree = M_PI*0.5 - colatitude;

                        i_lambda(lon_degree, lat_degree,
                                 grid_space_data[jlat * sphere2DDataConfig->grid_num_lon + ilon]);
                    }
                }

#endif
            }


            void DataGrid::grid_update_lambda_array(
                    std::function<void(int, int, double &)> i_lambda
            ) {

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD

#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS

                for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
                {
                    for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
                    {
                        i_lambda(i, j, grid_space_data[i*sphere2DDataConfig->grid_num_lat + j]);
                    }
                }

#else

                for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++) {
                    for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++) {
                        i_lambda(ilon, jlat, grid_space_data[jlat * sphere2DDataConfig->grid_num_lon + ilon]);
                    }
                }

#endif
            }


            void DataGrid::grid_update_lambda_array_idx(
                    std::function<void(int, double &)> i_lambda
            ) {
                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD

                for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++) {
                    i_lambda(i, grid_space_data[i]);
                }
            }


/*
 * Set values for all latitude and longitude degrees
 *
 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude sin(phi) \in [-1;1])
 */
            void DataGrid::grid_update_lambda_gaussian_grid(
                    std::function<void(double, double, double &)> i_lambda
            ) {
                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD

#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS

                for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
                {
                    double lon_degree = ((double)i/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

                    for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
                    {
                        double sin_phi = sphere2DDataConfig->lat_gaussian[j];

                        i_lambda(lon_degree, sin_phi, grid_space_data[i*sphere2DDataConfig->grid_num_lat + j]);
                    }
                }
#else

                for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++) {
                    double sin_phi = sphere2DDataConfig->lat_gaussian[jlat];

                    for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++) {
                        double lon_degree = ((double) ilon / (double) sphere2DDataConfig->grid_num_lon) * 2.0 * M_PI;

                        i_lambda(lon_degree, sin_phi, grid_space_data[jlat * sphere2DDataConfig->grid_num_lon + ilon]);
                    }
                }
#endif
            }


/*
 * Set values for all latitude and longitude degrees
 *
 * lambda function parameters:
 *   (longitude \in [0;2*pi], Cogaussian latitude cos(phi) \in [0;1])
 */
            void DataGrid::DataGrid::grid_update_lambda_cogaussian_grid(
                    std::function<void(double, double, double &)> i_lambda
            ) {

#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
                {
                    double lon_degree = (((double)i)/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

                    for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
                    {
                        double cos_phi = sphere2DDataConfig->lat_cogaussian[j];

                        /*
                         * IDENTITAL FORMULATION
                        double mu = shtns->ct[j];
                        double comu = sqrt(1.0-mu*mu);
                        */

                        i_lambda(lon_degree, cos_phi, grid_space_data[i*sphere2DDataConfig->grid_num_lat + j]);
                    }
                }
#else

                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++) {
                    double cos_phi = sphere2DDataConfig->lat_cogaussian[jlat];

                    for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++) {
                        double lon_degree = (((double) ilon) / (double) sphere2DDataConfig->grid_num_lon) * 2.0 * M_PI;

                        /*
                         * IDENTITAL FORMULATION
                        double mu = shtns->ct[j];
                        double comu = sqrt(1.0-mu*mu);
                        */

                        i_lambda(lon_degree, cos_phi, grid_space_data[jlat * sphere2DDataConfig->grid_num_lon + ilon]);
                    }
                }
#endif
            }


            void DataGrid::grid_update_lambda_sinphi_grid(
                    std::function<void(double, double, double &)> i_lambda
            ) {
                grid_update_lambda_gaussian_grid(i_lambda);
            }

            void DataGrid::grid_update_lambda_cosphi_grid(
                    std::function<void(double, double, double &)> i_lambda
            ) {
                grid_update_lambda_cogaussian_grid(i_lambda);
            }


/*
 * Set all values to zero
 */
            void DataGrid::grid_setZero() {
                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD

                for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
                    for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
                        grid_space_data[j * sphere2DDataConfig->grid_num_lon + i] = 0;
            }


/*
 * Set all values to a specific value
 */
            void DataGrid::grid_setValue(
                    double i_value
            ) {
                SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
                for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
                    for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
                        grid_space_data[j * sphere2DDataConfig->grid_num_lon + i] = i_value;
            }


/*
 * Set all values to a specific value
 */
            void DataGrid::grid_setValue(
                    int i_lon_idx,
                    int i_lat_idx,
                    double i_value
            ) {
#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS
                grid_space_data[i_lon_idx*sphere2DDataConfig->grid_num_lat + i_lat_idx] = i_value;
#else
                grid_space_data[i_lat_idx * sphere2DDataConfig->grid_num_lon + i_lon_idx] = i_value;
#endif
            }

            double DataGrid::grid_getValue(int i_lon_idx, int i_lat_idx) {
                return grid_space_data[i_lat_idx * sphere2DDataConfig->grid_num_lon + i_lon_idx];
            }


/**
 * Return the maximum error norm between this and the given data in physical space
 */
            double DataGrid::grid_reduce_max(
                    const DataGrid &i_sph_data
            ) {
                check(i_sph_data.sphere2DDataConfig);

                double error = -1;

                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++) {
                    error = std::max(
                            std::abs(
                                    grid_space_data[j] - i_sph_data.grid_space_data[j]
                            ),
                            error
                    );
                }
                return error;
            }


            double DataGrid::grid_reduce_rms() {
                double error = 0;

                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++) {
                    double &d = grid_space_data[j];
                    error += d * d;
                }

                return std::sqrt(error / (double) sphere2DDataConfig->grid_number_elements);
            }

            double DataGrid::grid_reduce_norm1() {
                double error = 0;

                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++) {
                    double &d = grid_space_data[j];
                    error += std::abs(d);
                }

                return error;
            }

            double DataGrid::grid_reduce_norm2() {
                double error = 0;

                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++) {
                    double &d = grid_space_data[j];
                    error += d * d;
                }

                return std::sqrt(error);
            }


            double DataGrid::grid_reduce_sum() const {
                double sum = 0;
                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++)
                    sum += grid_space_data[j];

                return sum;
            }


/**
 * return the maximum of all absolute values, use quad precision for reduction
 */
            double DataGrid::grid_reduce_sum_quad() const {
                double sum = 0;
                double c = 0;
#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum, c)
#endif
                for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++) {
                    double value = grid_space_data[i];

                    // Use Kahan summation
                    double y = value - c;
                    double t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }

                sum -= c;

                return sum;
            }


/**
 * Return the maximum error norm between this and the given data in physical space
 */
            double DataGrid::grid_reduce_max_abs(
                    const DataGrid &i_sph_data
            ) const {
                check(i_sph_data.sphere2DDataConfig);

                double error = -1;

                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++) {
                    error = std::max(
                            std::abs(
                                    grid_space_data[j] - i_sph_data.grid_space_data[j]
                            ),
                            error    // leave the error variable as the 2nd parameter. In case of NaN of the 1st parameter, std::max returns NaN
                    );
                }
                return error;
            }


/**
 * Return the maximum absolute value
 */
            double DataGrid::grid_reduce_max_abs() const {
                double error = -1;

                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++) {
                    error = std::max(
                            std::abs(grid_space_data[j]),
                            error        // leave the error variable as the 2nd parameter. In case of NaN of the 1st parameter, std::max returns NaN
                    );
                }
                return error;
            }


/**
 * Return the minimum value
 */
            double DataGrid::grid_reduce_min() const {
                double error = std::numeric_limits<double>::infinity();

                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++)
                    error = std::min(grid_space_data[j], error);

                return error;
            }


/**
 * Return the minimum value
 */
            double DataGrid::grid_reduce_max() const {
                double error = -std::numeric_limits<double>::infinity();

                for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++)
                    error = std::max(grid_space_data[j], error);

                return error;
            }


            bool DataGrid::grid_isAnyNaNorInf() const {
                for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++) {
                    if (std::isnan(grid_space_data[i]) || std::isinf(grid_space_data[i]) != 0)
                        return true;
                }

                return false;
            }


            void DataGrid::grid_print(
                    int i_precision
            ) const {
                if (i_precision >= 0)
                    std::cout << std::setprecision(i_precision);

                for (int j = (int) (sphere2DDataConfig->grid_num_lat - 1); j >= 0; j--) {
                    for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++) {
#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS
                        std::cout << grid_space_data[i*sphere2DDataConfig->grid_num_lat+j];
#else
                        std::cout << grid_space_data[j * sphere2DDataConfig->grid_num_lon + i];
#endif
                        if (i < sphere2DDataConfig->grid_num_lon - 1)
                            std::cout << "\t";
                    }
                    std::cout << std::endl;
                }
            }


            void DataGrid::grid_file_write(
                    const std::string &i_filename,
                    const char *i_title,
                    int i_precision
            ) const {
                std::ofstream file(i_filename, std::ios_base::trunc);

                if (i_precision >= 0)
                    file << std::setprecision(i_precision);

                file << "#TI " << i_title << std::endl;
                file << "#TX Longitude" << std::endl;
                file << "#TY Latitude" << std::endl;

                //file << "lat\\lon\t";
                // Use 0 to make it processable by python
                file << "0\t";

                for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++) {
//			double lon_degree = ((double)i/(double)sphere2DDataConfig->spat_num_lon)*2.0*M_PI;
                    double lon_degree = ((double) i / (double) sphere2DDataConfig->grid_num_lon) * 2.0 * M_PI;
                    lon_degree = lon_degree / M_PI * 180.0;

                    file << lon_degree;
                    if (i < sphere2DDataConfig->grid_num_lon - 1)
                        file << "\t";
                }
                file << std::endl;

                for (int j = sphere2DDataConfig->grid_num_lat - 1; j >= 0; j--) {
//        		double lat_degree =  M_PI*0.5 - acos(shtns->ct[j]);
                    double lat_degree = sphere2DDataConfig->lat[j];
                    lat_degree = lat_degree / M_PI * 180.0;

                    file << lat_degree << "\t";

                    for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++) {
#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS
                        file << grid_space_data[i*sphere2DDataConfig->grid_num_lat+j];
#else
                        file << grid_space_data[j * sphere2DDataConfig->grid_num_lon + i];
#endif
                        if (i < sphere2DDataConfig->grid_num_lon - 1)
                            file << "\t";
                    }
                    file << std::endl;
                }
                file.close();
            }


            void DataGrid::grid_file_write_lon_pi_shifted(
                    const char *i_filename,
                    const std::string &i_title,
                    int i_precision
            ) const {
                std::ofstream file(i_filename, std::ios_base::trunc);

                file << std::setprecision(i_precision);
                file << "#TI " << i_title << std::endl;
                file << "#TX Longitude" << std::endl;
                file << "#TY Latitude" << std::endl;

                //file << "lat\\lon\t";
                // Use 0 to make it processable by python
                file << "0\t";

                for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++) {
//			double lon_degree = ((double)i/(double)sphere2DDataConfig->spat_num_lon)*2.0*M_PI;
                    double lon_degree = ((double) i / (double) sphere2DDataConfig->grid_num_lon) * 2.0 * M_PI;
                    lon_degree = (lon_degree - M_PI) / M_PI * 180.0;

                    file << lon_degree;
                    if (i < sphere2DDataConfig->grid_num_lon - 1)
                        file << "\t";
                }
                file << std::endl;

                for (int j = sphere2DDataConfig->grid_num_lat - 1; j >= 0; j--) {
//        		double lat_degree =  M_PI*0.5 - acos(shtns->ct[j]);
                    double lat_degree = sphere2DDataConfig->lat[j];
                    lat_degree = lat_degree / M_PI * 180.0;

                    file << lat_degree << "\t";

                    for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++) {
                        int ia = i + sphere2DDataConfig->grid_num_lon / 2;
                        if (ia >= sphere2DDataConfig->grid_num_lon)
                            ia -= sphere2DDataConfig->grid_num_lon;

#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS
                        file << grid_space_data[ia*sphere2DDataConfig->grid_num_lat+j];
#else
                        file << grid_space_data[j * sphere2DDataConfig->grid_num_lon + ia];
#endif
                        if (i < sphere2DDataConfig->grid_num_lon - 1)
                            file << "\t";
                    }
                    file << std::endl;
                }
                file.close();
            }


            bool DataGrid::file_read_csv_grid(
                    const char *i_filename,
                    bool i_binary_data
            ) {

                ///Sphere2DDataConfig sphere2DDataConfig_ref;

                std::cout << "loading DATA from " << i_filename << std::endl;

                std::ifstream file(i_filename);

                for (int i = 0; i < 4; i++) {
                    std::string line;
                    std::getline(file, line);
                    std::istringstream iss(line);
                    std::vector<std::string> str_vector((std::istream_iterator<std::string>(iss)),
                                                        std::istream_iterator<std::string>());

                    if (i == 0) {
                        SWEET_ASSERT(str_vector.size() == 1);
                        SWEET_ASSERT(str_vector[0] == "#TI");
                    } else if (i == 1) {
                        SWEET_ASSERT(str_vector.size() == 2);
                        SWEET_ASSERT(str_vector[0] == "#TX");
                        SWEET_ASSERT(str_vector[1] == "Longitude");
                    } else if (i == 2) {
                        SWEET_ASSERT(str_vector.size() == 2);
                        SWEET_ASSERT(str_vector[0] == "#TY");
                        SWEET_ASSERT(str_vector[1] == "Latitude");
                    }
                        // Line 4: longitude values
                        // First character is "0"
                    else if (i == 3) {
                        SWEET_ASSERT((int) str_vector.size() == (int) sphere2DDataConfig->grid_num_lon + 1);
                        SWEET_ASSERT(str_vector[0] == "0");
                        for (int l = 0; l < sphere2DDataConfig->grid_num_lon; l++)
                                SWEET_ASSERT(std::abs(atof(str_vector[l + 1].c_str()) -
                                                      ((double) l / (double) sphere2DDataConfig->grid_num_lon) * 2.0 *
                                                      M_PI / M_PI * 180.0) < 1e-13);
                    }
                }


                ///sphere2DDataConfig_ref.setup(resx_ref, resy_ref, (int)((resx_ref * 2) / 3), (int)((resx_ref * 2) / 3), cart2DDataConfig->reuse_spectral_transformation_plans);
                ///*this = Sphere2DData_Grid(&sphere2DDataConfig_ref);


                for (int j = sphere2DDataConfig->grid_num_lat - 1; j >= 0; j--) {
                    std::string line;
                    std::getline(file, line);
                    std::istringstream iss(line);
                    std::vector<std::string> str_vector((std::istream_iterator<std::string>(iss)),
                                                        std::istream_iterator<std::string>());

                    if (!file.good()) {
                        std::cerr << "Failed to read data from file " << i_filename << " in line "
                                  << sphere2DDataConfig->grid_num_lat - 1 - j + 3 << std::endl;
                        return false;
                    }


                    SWEET_ASSERT((int) str_vector.size() == (int) sphere2DDataConfig->grid_num_lon + 1);
                    SWEET_ASSERT(
                            std::abs(atof(str_vector[0].c_str()) - sphere2DDataConfig->lat[j] / M_PI * 180.0) < 1e-13);

                    for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++) {
                        double val = atof(str_vector[i + 1].c_str());
#if SPHERE2D_DATA_GRID_LAYOUT == SPHERE2D_DATA_LAT_CONTIGUOUS
                        grid_space_data[i*sphere2DDataConfig->grid_num_lat+j] = val;
#else
                        grid_space_data[j * sphere2DDataConfig->grid_num_lon + i] = val;
#endif
                    }
                }

                file.close();
                std::cout << "DATA loaded OK" << std::endl;

                return true;
            }

            void DataGrid::file_write_raw(
                    const std::string &i_filename
            ) const {
                std::fstream file(i_filename, std::ios::out | std::ios::binary);
                file.write((const char *) grid_space_data, sizeof(double) * sphere2DDataConfig->grid_number_elements);
            }


            void DataGrid::file_read_raw(
                    const std::string &i_filename
            ) const {
                std::fstream file(i_filename, std::ios::in | std::ios::binary);
                file.read((char *) grid_space_data, sizeof(double) * sphere2DDataConfig->grid_number_elements);
            }


            void DataGrid::print_debug(
                    const char *name
            ) const {
                std::cout << name << ":" << std::endl;
                std::cout << "                min: " << grid_reduce_min() << std::endl;
                std::cout << "                max: " << grid_reduce_max() << std::endl;
                std::cout << "                sum: " << grid_reduce_sum() << std::endl;
                std::cout << std::endl;
            }


            void DataGrid::print() const {
                for (std::size_t idx = 0; idx < sphere2DDataConfig->grid_number_elements; idx++)
                    std::cout << grid_space_data[idx] << "\t";
                std::cout << std::endl;
            }


        }
    }
}
