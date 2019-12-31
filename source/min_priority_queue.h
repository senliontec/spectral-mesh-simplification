#pragma once

#include <vector>
#include <functional>

template<typename T, typename IndexFn, typename CompareFn = std::less<T>>
class min_priority_queue
{
public:
	min_priority_queue(min_priority_queue&&) = default;
	min_priority_queue(const min_priority_queue&) = default;
	min_priority_queue& operator=(min_priority_queue&&) = default;
	min_priority_queue& operator=(const min_priority_queue&) = default;
	min_priority_queue(IndexFn idx = IndexFn(), CompareFn cmp = CompareFn()) :
		compare_less(cmp),
		index_of(idx)
	{
	}

	/// Return the number of elements
	size_t size() const;

	/// Return whether there are no elements
	[[nodiscard]] bool empty() const;

	/// Clear both elements and indices, without deallocating the memory
	void clear();

	/// Reserve memory, with a different capacity for indices (usually bigger)
	void reserve(size_t num_elements, size_t num_indices);

	/// Reserve memory, same capacity for elements and indices
	void reserve(size_t num_elements);

	/// Return the smallest element w.r.t. the comparison predicate
	const T& top() const;

	/// Add an element to the queue, or update it if already present
	void push(const T& data);

	/// Add an element to the queue, or update it if already present (move version)
	void push(T&& data);

	/// Pop the first element (which is the smallest one) by moving it in the returned value
	T pop();

	/// Update the element at the specified index, if it exists (return whether it was updated)
	template<typename F>
	bool update(size_t index, F updater);

	/// Return whether there is an element associated with the specified index in the queue
	bool contains(size_t idx) const;

private:
	void ensure_index(size_t idx);
	size_t right_child(size_t i) const;
	size_t left_child(size_t i) const;
	size_t parent(size_t i) const;
	void swap(size_t e1, size_t e2);
	void push_down(size_t e);
	void push_up(size_t e);

private:
	std::vector<T> elements;
	std::vector<size_t> indices;
	CompareFn compare_less;
	IndexFn index_of;
};

template<typename T, typename IndexFn, typename CompareFn>
size_t min_priority_queue<T, IndexFn, CompareFn>::size() const
{
	return elements.size();
}

template<typename T, typename IndexFn, typename CompareFn>
bool min_priority_queue<T, IndexFn, CompareFn>::empty() const
{
	return elements.empty();
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::clear()
{
	elements.clear();
	indices.clear();
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::reserve(size_t num_elements, size_t num_indices)
{
	elements.reserve(num_elements);
	indices.resize(num_indices, size_t(-1));
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::reserve(size_t num_elements)
{
	reserve(num_elements, num_elements);
}

template<typename T, typename IndexFn, typename CompareFn>
const T& min_priority_queue<T, IndexFn, CompareFn>::top() const
{
	return elements.front();
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::push(const T& data)
{
	const size_t idx = index_of(data);
	ensure_index(idx);
	size_t e = indices[idx];
	if(e != size_t(-1))
	{
		elements[e] = data;
		push_down(e);
	}
	else
	{
		e = elements.size();
		indices[idx] = e;
		elements.push_back(data);
	}
	push_up(e);
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::push(T&& data)
{
	const size_t idx = index_of(data);
	ensure_index(idx);
	size_t e = indices[idx];
	if(e != size_t(-1))
	{
		elements[e] = std::move(data);
		push_down(e);
	}
	else
	{
		e = elements.size();
		indices[idx] = e;
		elements.push_back(std::move(data));
	}
	push_up(e);
}

template<typename T, typename IndexFn, typename CompareFn>
T min_priority_queue<T, IndexFn, CompareFn>::pop()
{
	T p = std::move(elements[0]);
	swap(0, elements.size() - 1);
	indices[index_of(p)] = size_t(-1);
	elements.pop_back();
	push_down(0);
	return p;
}

template<typename T, typename IndexFn, typename CompareFn>
template<typename F>
bool min_priority_queue<T, IndexFn, CompareFn>::update(size_t index, F op)
{
	const size_t el = index < indices.size() ? indices[index] : size_t(-1);
	if(el == size_t(-1)) return false;
	op(elements[el]);
	push_down(el);
	push_up(el);
	return true;
}

template<typename T, typename IndexFn, typename CompareFn>
bool min_priority_queue<T, IndexFn, CompareFn>::contains(size_t idx) const
{
	return idx < indices.size() && indices[idx] != size_t(-1);
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::ensure_index(size_t idx)
{
	if(idx < indices.size()) return;
	size_t s = std::max(size_t(16), indices.size());
	/// TODO: assert(idx < size_t(-1) / 2)
	while(s <= idx) s *= 2;
	indices.resize(s, size_t(-1));
}

template<typename T, typename IndexFn, typename CompareFn>
size_t min_priority_queue<T, IndexFn, CompareFn>::right_child(size_t i) const
{
	return 2 * i + 2;
}

template<typename T, typename IndexFn, typename CompareFn>
size_t min_priority_queue<T, IndexFn, CompareFn>::left_child(size_t i) const
{
	return 2 * i + 1;
}

template<typename T, typename IndexFn, typename CompareFn>
size_t min_priority_queue<T, IndexFn, CompareFn>::parent(size_t i) const
{
	return (i - 1) / 2;
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::swap(size_t e1, size_t e2)
{
	std::swap(indices[index_of(elements[e1])], indices[index_of(elements[e2])]);
	std::swap(elements[e1], elements[e2]);
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::push_down(size_t e)
{
	bool pushed = false;
	do
	{
		size_t min = e;
		const size_t l = left_child(e), r = right_child(e);
		if(l < elements.size() && compare_less(elements[l], elements[min])) min = l;
		if(r < elements.size() && compare_less(elements[r], elements[min])) min = r;
		pushed = min != e;
		if(pushed)
		{
			swap(e, min);
			e = min;
		}
	}
	while(pushed);
}

template<typename T, typename IndexFn, typename CompareFn>
void min_priority_queue<T, IndexFn, CompareFn>::push_up(size_t e)
{
	while(e > 0 && compare_less(elements[e], elements[parent(e)]))
	{
		swap(e, parent(e));
		e = parent(e);
	}
}
