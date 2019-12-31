#pragma once

#include "mesh.h"
#include <charconv>
#include "util.h"

struct reduction_metric
{
	using serialized_parameters = std::vector<std::pair<std::string, std::string>>;

	virtual ~reduction_metric()
	{
	}

	virtual serialized_parameters save() const
	{
		return {};
	}

	virtual void load(const serialized_parameters&)
	{
	}

	virtual void setup(const mesh&, const half_edge_connectivity&)
	{
	}

	virtual bool collapse_still_valid(uint32_t /* half_edge */) const
	{
		return true;
	}

	virtual bool flip_still_valid(uint32_t /* half_edge */) const
	{
		return true;
	}

	virtual std::pair<Eigen::Vector3d, double> cost_collapse(uint32_t /* half_edge */, uint32_t /* to_keep */, uint32_t /* to_remove */) const
	{
		return { Eigen::Vector3d::Zero(), std::numeric_limits<double>::quiet_NaN() };
	}

	virtual unsigned int post_collapse(uint32_t /* half_edge */, uint32_t /* to_keep */, uint32_t /* to_remove */)
	{
		return 0;
	}

	virtual double cost_flip(uint32_t /* half_edge */) const
	{
		return std::numeric_limits<double>::quiet_NaN();
	}

	virtual unsigned int post_flip(uint32_t /* half_edge */)
	{
		return 0;
	}

protected:
	template<typename T>
	static void serialize(serialized_parameters& p, const char* name, const T& value)
	{
		if constexpr(std::is_same_v<T, bool>)
			p.emplace_back(name, value ? "true" : "false");
		else
			p.emplace_back(name, std::to_string(value));
	}

	template<typename T>
	static void deserialize(const serialized_parameters& p, const char* name, T& value)
	{
		if constexpr(std::is_same_v<T, bool>)
		{
			for(const std::pair<std::string, std::string>& x : p)
				if(x.first == name)
					value = ci_equal(x.second.c_str(), "true")
						|| ci_equal(x.second.c_str(), "yes")
						|| ci_equal(x.second.c_str(), "1")
						|| ci_equal(x.second.c_str(), "on");
		}
		else
		{
			for(const std::pair<std::string, std::string>& x : p)
				if(x.first == name)
					std::from_chars(x.second.data(), x.second.data() + x.second.size(), value);
		}
	}
};
