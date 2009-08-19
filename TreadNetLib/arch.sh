#!/bin/sh

system=`uname -s``uname -r`

case $system in
        SunOS4*)
                system=sunOS
                ;;
        SunOS5*)
                system=solaris
                ;;
        IRIX*)
                system=irix
		;;
        HP-UX*)
                system=hpux
                ;;
        Linux*)
                system=linux
                ;;
        OSF1*)
                system=osf1
                ;;
        *NT*)
                system=nt
                ;;
        BSD*)
                system=BSD
                ;;
        NetBSD*)
                system=netbsd
                ;;
        FreeBSD*)
                system=FreeBSD
                ;;

        *)
                echo "Unknown system \"$system\"."
                exit 1
                ;;              
esac

machine=`uname -m`
case $machine in
        sun4*)
                machine=sparc
                ;;
        sparc*)
                machine=sparc
                ;;
        9000/*)
                machine=hp
                ;;
        i[3456]86 | i86pc | x86_64)
                machine=intel
                ;;
        alpha)
                machine=dec
                ;;
        IP*)
		processor=`hinv -t cpu`
		case $processor in
		    *R12000*)
			machine=mips4
			;;
		    *R10000*)
			machine=mips4
			;;
		    *R4400*)
		        machine=mips3
			;;
		esac
                ;;
        *)
                echo "Unknown machine \"$machine\"."
                exit 1
                ;;              
esac    
 
if test x`hostname` = xmerlin
then
        machine=hp
        system=pci
fi

if test x`hostname` = xpluto
then
        machine=hp
        system=pci
fi

if test x`hostname` = xkitty
then
        machine=hp
        system=dino
fi

 
echo "$machine"_"$system"
exit 0
