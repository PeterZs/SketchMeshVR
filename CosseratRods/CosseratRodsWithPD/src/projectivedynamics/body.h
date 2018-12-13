#ifndef __BODY_PD__
#define __BODY_PD__

class Body {
public:
	Body() {}
	virtual void Integrate() {}
	virtual void Render() {}
private:
};

#endif //__BODY_PD__