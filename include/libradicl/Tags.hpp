
#ifndef LIBRADICL_TAGS_HPP
#define LIBRADICL_TAGS_HPP



#include "Type.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>


namespace RAD
{


class Buffer;


class Tag_Defn
{
private:

    typedef std::pair<Type::str, Type::u8> tag_desc;

    std::vector<tag_desc> file_tag;
    std::vector<tag_desc> read_tag;
    std::vector<tag_desc> aln_tag;


public:

    template <typename T_>
    void add_file_tag(const std::string& name) { static_assert(is_RAD_type<T_>()); file_tag.emplace_back(name, T_::type_id()); }

    template <typename T_>
    void add_read_tag(const std::string& name) { static_assert(is_RAD_type<T_>()); read_tag.emplace_back(name, T_::type_id()); }

    template <typename T_>
    void add_aln_tag(const std::string& name)  { static_assert(is_RAD_type<T_>()); aln_tag.emplace_back(name, T_::type_id()); }

    void write(Buffer& buf) const;
};


class Tag_List
{
private:

    typedef Type::variant_t Tag;
    std::vector<Tag> tag;


    template <typename T_sink_>
    class Writer : public boost::static_visitor<>
    {
    private:

        T_sink_& sink;

    public:

        explicit Writer(T_sink_& sink): sink(sink) {}

        template <typename T_>
        void operator()(const T_& operand) const { sink.add(operand); }
    };


public:

    void add(const Tag& t) { tag.emplace_back(t); }

    void clear() { tag.clear(); }

    template <typename T_sink_>
    void write(T_sink_& sink) const;
};


template <typename T_sink_>
inline void Tag_List::write(T_sink_& sink) const
{
    Writer<T_sink_> writer(sink);
    std::for_each(tag.cbegin(), tag.cend(), boost::apply_visitor(writer));
}

}



#endif
