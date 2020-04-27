
	function dcosd(arg)
	implicit none
	double precision arg, dcosd

	dcosd = dcos(arg/1.8D02*3.141592653589793238462643383279D0)
c	dcosd = dcos(arg/1.8D02*acos(-1.0D0))
	return
	end

	function dsind(arg)
	implicit none
	double precision arg, dsind

	dsind = dsin(arg/1.8D02*3.141592653589793238462643383279D0)
c	dsind = dsin(arg/1.8D02*acos(-1.0D0))
	return
	end

	function dacosd(arg)
	implicit none
	double precision arg, dacosd

	dacosd = dacos(arg)/3.141592653589793238462643383279D0*1.8D02
c	dacosd = dacos(arg)/acos(-1.0D0)*1.8D02
	return
	end

	function dasind(arg)
	implicit none
	double precision arg, dasind

	dasind = dasin(arg)/3.141592653589793238462643383279D0*1.8D02
c	dasind = dasin(arg)/acos(-1.0D0)*1.8D02
	return
	end

	function dtand(arg)
	implicit none
	double precision arg, dtand

	dtand = dtan(arg/1.8D02*3.141592653589793238462643383279D0)
c	dtand = dtan(arg/1.8D02*acos(-1.0D0))
	return
	end

	function datand(arg)
	implicit none
	double precision arg, datand

	datand = datan(arg)/3.141592653589793238462643383279D0*1.8D02
c	datand = datan(arg)/acos(-1.0D0)*1.8D02
	return
	end
