#include "stdafx.h"
#include "XMLDoc.h"

#include <assert.h>

using namespace std;
//using namespace boost;

#include "tinyxml\tinyxml.h"

namespace tixml{

	string trim_copy(string str){
		int start = 0, end = str.size() - 1;
		while(str[start] == ' ' || str[start] == '\n' || str[start] == '\t' || str[start] == '\r') start ++;
		while(str[end] == ' ' || str[end] == '\n' || str[end] == '\t' || str[end] == '\r') end --;
		return str.substr(start, end-start+1);
	}

	void split(vector<string>& keys, string key, char splitor){
		keys.clear();
		int begin = 0;
		for(int i = 0; i < (int)key.size(); i++){
			if(key[i] == splitor){
				if(i - 1 < begin) keys.push_back("");
				else keys.push_back(key.substr(begin, i - begin));
				begin = i + 1;
				i ++;
			}
		}//cout<<begin<<" "<<key.size();
		keys.push_back(key.substr(begin, key.size() - begin));
	}

	static void convertFromTiXML(XMLNode&node,TiXmlElement* elem){
		if(elem==NULL)
			return;

		node.setName(elem->Value());		
		if(elem->GetText()){
			node.setText(string(elem->GetText()));
		} else 
			node.setText(string(""));

		TiXmlAttribute* tiAttrib = elem->FirstAttribute();
		while(tiAttrib){
			node.attribs()[tiAttrib->Name()] = tiAttrib->Value();
			tiAttrib = tiAttrib->Next();
		}

		TiXmlElement* tiChild = elem->FirstChildElement();
		while(tiChild){

			TiXmlElement* tiChildElem = (tiChild);
			if(  tiChildElem ){
				XMLNode* child= new XMLNode();
				convertFromTiXML(*child,tiChildElem);
				node.children().push_back(child);
			}

			tiChild = tiChild->NextSiblingElement();
		}
	}

	static void convertToTiXML(XMLNode* pNode,TiXmlElement* elem){
		if(pNode==NULL)
			return;

		elem->SetValue(pNode->name());
		for(std::map<std::string,std::string>::iterator it = pNode->attribs().begin(); it!=pNode->attribs().end();it++){
			elem->SetAttribute(it->first,it->second);
		}

		string trimmedText = trim_copy(pNode->text<string>());
		if(trimmedText.length()>0) {
			TiXmlText * tiText = new TiXmlText(trimmedText);
			elem->LinkEndChild(tiText);
		}

		for(unsigned int i=0; i<pNode->children().size();i++){
			TiXmlElement* childElem = new TiXmlElement("");
			convertToTiXML(pNode->children()[i],childElem);
			elem->LinkEndChild(childElem);
		}

	}

	XMLDoc::XMLDoc(void)
	{
		m_pRoot = NULL;
	}


	XMLDoc::~XMLDoc(void)
	{
		if(m_pRoot)
			delete m_pRoot;
	}

	void XMLDoc::Create(const std::string & rootName,const std::map<std::string,std::string> & rootAttribs) {
		if(m_pRoot)
			delete m_pRoot;
		m_pRoot = new XMLNode(rootName,"",rootAttribs);
	}
	void XMLDoc::Create(const std::string & rootName){
		std::map<std::string,std::string> rootAttribs;
		Create(rootName,rootAttribs);
	}

	bool XMLDoc::Load(const char* fn) {
		
		TiXmlDocument doc;
		if(!doc.LoadFile(fn))
			return false;
		if(m_pRoot)
			delete m_pRoot;
		m_pRoot = new XMLNode();
		TiXmlElement* rootElem = doc.RootElement();
		convertFromTiXML(*m_pRoot,rootElem);

		return true;
	}
	bool XMLDoc::Save(const char* fn) {

		TiXmlDocument doc;
		TiXmlDeclaration * tiDecl = new TiXmlDeclaration("1.0","","");
		TiXmlElement * tiRoot = new TiXmlElement("");
		convertToTiXML(getRoot(),tiRoot);

		doc.LinkEndChild(tiDecl);
		doc.LinkEndChild(tiRoot);

		return doc.SaveFile(fn);

		return false;
	}

	void XMLDoc::test(){

		//load
		XMLDoc doc;
		if(!doc.Load("test.xml")){
			cout<<"'test.xml' is missing. exit!\n";
			return;
		}
		cout<<"id.@name = "<<doc.get<string>("id.@name", 3)<<endl;
		cout<<"id4 = "<<doc.get<string>("id", 4)<<endl;

		//save
		tixml::XMLDoc doc2;
		doc2.Create("root");
		tixml::XMLNode* pNode1 = new tixml::XMLNode("name", "123");
		doc2.getRoot()->children().push_back(pNode1);
		for(int i = 0; i < 5; i++){
			tixml::XMLNode * pChar = new tixml::XMLNode("id", "hi");
			pChar->setAttrib<int>("name", i);
			pNode1->children().push_back(pChar);
		}

		doc2.Save("test2.xml");
	}
	
	unsigned int XMLNodeList::length() const {
		return m_list.size();
	}



	XMLNodeList XMLNodeList::getChildNodeList(const std::string& name) const{
		vector<XMLNode* > childList;
		for(unsigned int i=0;i<length();i++){
			std::vector<XMLNode*>& list = m_list[i]->children();
			for(unsigned int j=0;j<list.size();j++){
				if(list[j]->name() == name){
					childList.push_back(list[j]);
				}
			}			
		}
		return XMLNodeList(childList);
	}

	XMLNodeList XMLNodeList::get(const std::string& key) const {
		
		vector<string> keys;
		split(keys, key, '.');
		//split(keys,trim_copy(key),is_any_of("."));
		
		XMLNodeList res(*this);
		for(unsigned int i=0;i<keys.size();i++) {
			res = res.getChildNodeList(keys[i]);
		}		
		return res;
	}
}
