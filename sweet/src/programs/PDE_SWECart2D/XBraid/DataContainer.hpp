#ifndef PROGRAMS_PDE_SWECART2D_XRAID_DATACONTAINER_HPP
#define PROGRAMS_PDE_SWECART2D_XRAID_DATACONTAINER_HPP

#include <sweet/XBraid/Vector.hpp>
#include <programs/PDE_SWECart2D/DataContainer/Simulation.hpp>


namespace PDE_SWECart2D {
namespace XBraid {

class DataContainer :
	public sweet::XBraid::Vector
{

public:
	sweet::Data::Cart2D::Config*				config = nullptr;
	PDE_SWECart2D::DataContainer::Simulation*		data = nullptr;

public:
	DataContainer()
	{
	}

public:
	~DataContainer()
	{
		clear();
	}

public:
	DataContainer(
		const DataContainer &i_vector
	)
	{
		config = i_vector.config;
		data->op_setVector(*i_vector.data);
		level = i_vector.level;
		N = i_vector.N;
		var_names = i_vector.var_names;
	};

public:
	DataContainer& operator=(
		const DataContainer &i_vector
	)
	{
		config = i_vector.config;
		data->op_setVector(*i_vector.data);
		level = i_vector.level;
		N = i_vector.N;
		var_names = i_vector.var_names;

		return *this;
	}


public:
	void clear()
	{
		if (data)
		{
			data->clear();
			data = nullptr;
		}
	}

public:
	void setup(
			sweet::Data::Cart2D::Config* i_config,
			int i_level
	)
	{
		config = i_config;
		level = i_level;
		allocate_data();
		N = data->N;
		var_names = data->var_names;
	}

private:
	static inline
	DataContainer* cast(sweet::XBraid::Vector* i_U)
	{
		return static_cast<DataContainer*>(i_U);
	}

private:
	static inline
	DataContainer& cast(sweet::XBraid::Vector &i_U)
	{
		return static_cast<DataContainer&>(i_U);
	}

private:
	static inline
	const DataContainer& cast(const sweet::XBraid::Vector &i_U)
	{
		return static_cast<const DataContainer&>(i_U);
	}


public:
	void op_setVector(
			const sweet::XBraid::Vector& i_U
	)	override
	{
		const DataContainer &U = cast(i_U);
		data->op_setVector(*U.data);
	}

public:
	void op_addVector(
			const sweet::XBraid::Vector &i_U
	)	override
	{
		const DataContainer &U = cast(i_U);
		data->op_addVector(*U.data);
	}

public:
	void op_subVector(
			const sweet::XBraid::Vector &i_U
	)	override
	{
		const DataContainer &U = cast(i_U);
		data->op_subVector(*U.data);
	}


public:
	void op_mulScalar(
			const double i_value
	)	override
	{
		data->op_mulScalar(i_value);
	}

public:
	virtual
	void allocate_data() override
	{
		data = new PDE_SWECart2D::DataContainer::Simulation;

		data->setup(this->config);
	}

	void restrict(sweet::XBraid::Vector* i_data) override
	{

		PDE_SWECart2D::XBraid::DataContainer* data2 = cast(i_data);
		data->restrict(*data2->data);
	}

	void pad_zeros(sweet::XBraid::Vector* i_data) override
	{

		PDE_SWECart2D::XBraid::DataContainer* data2 = cast(i_data);
		data->pad_zeros(*data2->data);
	}

	void op_setVector(sweet::XBraid::Vector* i_data) override
	{

		PDE_SWECart2D::XBraid::DataContainer* data2 = cast(i_data);
		data->op_setVector(*data2->data);
	}


	double reduceSum() const override
	{
		return data->reduceSum();
	}

	double reduceMaxAbs() const override
	{
		return data->reduceMaxAbs();
	}

	double reduceSum(int i) const override
	{
		return data->reduceSum(i);
	}

	double reduceMaxAbs(int i) const override
	{
		return data->reduceMaxAbs(i);
	}

	double reduceMaxAbs(int i, int rnorm) const override
	{
		return data->reduceMaxAbs(i, rnorm);
	}

	double reduceNormL1Grid(bool normalized = false) const override
	{
		return data->reduceNormL1Grid(normalized);
	}

	double reduceNormL2Grid(bool normalized = false) const override
	{
		return data->reduceNormL2Grid(normalized);
	}

	double reduceNormLinfGrid() const override
	{
		return data->reduceNormLinfGrid();
	}

	double reduceNormL1Grid(int i, bool normalized = false) const override
	{
		return data->reduceNormL1Grid(i, normalized);
	}

	double reduceNormL2Grid(int i, bool normalized = false) const override
	{
		return data->reduceNormL2Grid(i, normalized);
	}

	double reduceNormLinfGrid(int i) const override
	{
		return data->reduceNormLinfGrid(i);
	}

	void serialize(std::complex<double> *i_data) override
	{
		return data->serialize(i_data);
	}

	void deserialize(std::complex<double> *i_data) override
	{
		return data->deserialize(i_data);
	}

	void fileLoad(
			const char* i_filename_template,
			std::string i_path,
			std::string i_output_file_mode,
			double i_t,
			double i_timescale
	) override
	{
		for (int ivar = 0; ivar < data->N; ivar++)
		{
			std::string var_name = data->var_names[ivar];
			char buffer[1024];
			sprintf(buffer, i_filename_template, var_name.c_str(), i_t * i_timescale);
			std::string buffer2 = i_path + "/" + std::string(buffer);
			data->fileLoad(buffer2.c_str(), i_output_file_mode);
		}
	}
};

}}

#endif
