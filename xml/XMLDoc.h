/********************************************************************
	Copyright:		TNList, Tsinghua Univ.
	Author:			Li-Qian Ma
	Date:			2013-05-13
	Description:	A tiny xml parser.
********************************************************************/

#pragma once

#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include <sstream>
//#include <boost/shared_ptr.hpp>
//#include <helper/helper.h>

namespace tixml{
	
	template<typename T> std::string toStr(const T& t){
		std::stringstream ss;
		ss << t;
		return ss.str();
	}

	template<typename T> T fromStr(const std::string&str){
		T t;
		std::stringstream ss(str);
		ss >> t;
		return t;
	}


	class XMLNode;
	class XMLNodeList;

	class XMLNodeList{
	public:

		unsigned int length() const;

		XMLNodeList get(const std::string& key) const;

		XMLNode& operator [] (unsigned int i) const{
			assert(i>=0 && i<length());
			return *(m_list[i]);
		}

		XMLNode& operator () (const std::string& key, unsigned int i) const{
			return get(key)[i];
		}

		template <typename T> std::vector<T> get(const std::string& keyAndAttrib) const{			

			std::string::size_type pos = keyAndAttrib.find_last_of("@");
			
			if(pos!=std::string::npos){
				std::string key = keyAndAttrib.substr(0,pos-1);
				std::string attrib = keyAndAttrib.substr(pos+1);

				XMLNodeList list = get(key);

				std::vector<T> res(list.length());
				for(unsigned int i=0;i<list.length();i++){
					res[i] = list.getNode(i)->attrib<T>(attrib);
				}
				return res;
			} else {
				XMLNodeList list = get(keyAndAttrib);
				std::vector<T> res(list.length());
				for(unsigned int i=0;i<list.length();i++){
					res[i] = list.getNode(i)->text<T>();
				}
				return res;
			}
		}
		template <typename T> T get(const std::string& keyAndAtrrib,unsigned int i) const{			
			return get<T>(keyAndAtrrib)[i];
		}
		template <typename T> T unique(const std::string& keyAndAttrib) const{			
			return get<T>(keyAndAttrib,0);
		}


	public: 
		XMLNodeList();
		XMLNodeList(const XMLNodeList& list){
			m_list = list.m_list;
		}
		XMLNodeList(const std::vector<XMLNode* > & list){
			m_list = list;
		}

		void removeNode(unsigned int i) {
			m_list.erase(m_list.begin() + i);
		}
		XMLNode* getNode(unsigned int i) {
			assert(i>=0 && i<length());
			return m_list[i];
		}
		void addNode(XMLNode* newNode){
			m_list.push_back(newNode);
		}

	protected:
		XMLNodeList getChildNodeList(const std::string& name) const;

		std::vector<XMLNode* > m_list;

	};

	class XMLNode{
	public:

		XMLNode(){

		}
		XMLNode(const std::string& name, const std::string& text){
			m_name = name;
			m_text = text;
		}
		XMLNode(const std::string& name, const std::string& text, const std::map<std::string,std::string> & attribs){
			m_name = name;
			m_text = text;
			m_attribs = attribs;
		}

		~XMLNode(){
			for(unsigned int i=0;i<m_children.size();i++){
				if(m_children[i])
					delete m_children[i];
			}
		}

		template<typename T> T attrib(const std::string& name) const {
			return fromStr<T>(m_attribs.find(name)->second);
		}
		template<typename T> T text() const {
			return fromStr<T>(m_text);
		}
		const std::string& name() const {
			return m_name;
		}
		template<typename T> void setAttrib(const std::string& name,const T& val){
			m_attribs[name] = toStr(val);
		}
		template<typename T> void setText(const T& val){
			m_text = toStr(val);
		}
		void setName(const std::string& name) {
			m_name = name;
		}

		const std::vector<XMLNode*>& children() const {return m_children;}
		std::vector<XMLNode*>& children(){ return m_children;}

		const std::map<std::string,std::string> & attribs() const {
			return m_attribs;
		}
		std::map<std::string,std::string> & attribs() {
			return m_attribs;
		}

		XMLNodeList asNodeList(){
			std::vector<XMLNode* > list;
			list.push_back(this);
			return XMLNodeList(list);
		}

	public:
		XMLNodeList getRecurNodeList(std::string& key) const;

	protected:
		std::vector<XMLNode*> m_children;
		std::string m_name;		
		std::string m_text;
		std::map<std::string,std::string> m_attribs;
	};

	class XMLDoc
	{
	public:
		XMLDoc(void);
		~XMLDoc(void);

		static void test();

		bool Load(const char* fn);
		bool Save(const char* fn);
		void Create(const std::string & rootName,const std::map<std::string,std::string> & rootAttribs);
		void Create(const std::string & rootName);

		XMLNodeList get(const std::string& key) const {
			return m_pRoot->asNodeList().get(key);
		}

		XMLNode& operator () (const std::string& key, unsigned int i) const{
			return get(key)[i];
		}

		template <typename T> std::vector<T> get(const std::string& keyAndAtrrib) const{
			return m_pRoot->asNodeList().get<T>(keyAndAtrrib);
		}
		template <typename T> T get(const std::string& keyAndAttrib,unsigned int i) const{			
			return m_pRoot->asNodeList().get<T>(keyAndAttrib,i);
		}
		template <typename T> T unique(const std::string& keyAndAttrib) const{			
			return get<T>(keyAndAttrib,0);
		}

		XMLNode* getRoot() {
			return m_pRoot;
		}
	protected:
		XMLNode* m_pRoot;		
	};


}